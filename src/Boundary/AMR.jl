#=
AMR for boundaries, especially for immersed boundaries.
=#

function solid_static_refine_flag(forest, which_tree, quadrant)::Cint
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        ka = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ibs = ka.kinfo.config.IB
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(P4estPsData, qp.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        ps_data.bound_enc<0 && return Cint(1)
        for i in eachindex(ibs)
            ib = ibs[i]
            if search_radius_flag(ib, ps_data.midpoint, ps_data.ds)
                ps_data.bound_enc = i
                return Cint(1)
            end
        end
        return Cint(0)
    end
end
function solid_static_coarsen_flag(forest::Ptr{p8est_t}, which_tree, quadrants)::Cint
    GC.@preserve forest which_tree quadrants begin
        DIM = 3
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p8est_quadrant_t}}, quadrants, 2^DIM)
        for i = 1:(2^DIM)
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            !isa(ps_data, InsideSolidData)&&return Cint(0)
        end
        return Cint(1)
    end
end
function solid_static_replace!(
    forest,
    which_tree,
    num_out,
    out_quads::Ptr{T},
    num_in,
    in_quads,
)::Cvoid where {T}
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        ka = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{T}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{T}, in_quads, num_in)
        solid_static_replace!(
            Val(Int(num_out)),
            out_quads_wrap,
            in_quads_wrap,
            which_tree,
            ka,
        )
        return nothing
    end
end
function solid_static_replace!(
    ::Val{1},
    out_quad,
    in_quads,
    which_tree,
    ka::KA{DIM,NDF},
) where {DIM,NDF}# refine replace
    fp = PointerWrapper(ka.kinfo.forest.p4est)
    trees = ka.kdata.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_out_quad = PointerWrapper(out_quad[1])
    Odata = unsafe_pointer_to_objref(
        pointer(PointerWrapper(P4estPsData, pw_out_quad.p.user_data[]).ps_data),
    )
    index = findfirst(x -> x === Odata, datas)
    deleteat!(datas, index)
    if isa(Odata, InsideSolidData) #
        ib = ka.kinfo.config.IB[-Odata.bound_enc]
        for i = 1:(2^DIM)
            pw_in_quad = PointerWrapper(in_quads[i])
            dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
            # ds = Odata.ds*0.5
            # midpoint = Odata.midpoint+ 0.5 * ds .* RMT[DIM][i]
            ds, midpoint = quad_to_cell(fp, which_tree, pw_in_quad)
            if ghost_cell_flag(ib, midpoint, ds) # New born ghost cell
                ps_data = PsData(DIM, NDF; bound_enc = Odata.bound_enc, midpoint, ds)
                ps_data.vs_data = initialize_solid_vs_data(ib, ka.kinfo)
            else
                ps_data = InsideSolidData{DIM,NDF}(Odata.bound_enc, midpoint, ds)
            end
            insert!(datas, index - 1 + i, ps_data)
            dp[] = P4estPsData(pointer_from_objref(ps_data))
        end
    elseif Odata.bound_enc<0
        ib = ka.kinfo.config.IB[-Odata.bound_enc]
        flag = Odata.neighbor.state[1]==BALANCE_FLAG
        for i = 1:(2^DIM)
            pw_in_quad = PointerWrapper(in_quads[i])
            dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
            if flag # Balance flag
                ps_data = Odata.neighbor.data[1][i]
            else  # Not-a-balance
                ds = Odata.ds*0.5
                midpoint = @. Odata.midpoint+0.5 * ds * RMT[DIM][i]
                if solid_flag(ib, midpoint) # Inside solid
                    if ghost_cell_flag(ib, midpoint, ds)
                        ps_data = ps_copy(Odata)
                        ps_data.bound_enc = Odata.bound_enc
                        ps_data.ds .*= 0.5
                        vs_data = ps_data.vs_data
                        @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
                        if i==1
                            ps_data.neighbor.state[2] = BALANCE_FLAG
                            ps_data.neighbor.data[2] = [Odata]
                        end
                    else
                        ps_data = InsideSolidData{DIM,NDF}(Odata.bound_enc, midpoint, ds)
                    end
                else # fluid cell
                    ps_data = ps_copy(Odata)
                    ps_data.bound_enc = -Odata.bound_enc
                    ps_data.ds .*= 0.5
                    vs_data = ps_data.vs_data
                    @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
                    for k in axes(vs_data.df, 2)
                        for j in axes(vs_data.df, 1)
                            @views vs_data.df[j, k]+=min(
                                abs(
                                    vs_data.df[j, k]/(
                                        0.5*dot(ps_data.ds, abs.(vs_data.sdf[j, k, :]))+EPS
                                    ),
                                ),
                                1.0,
                            )*dot(vs_data.sdf[j, k, :], 0.5 * ps_data.ds .* RMT[DIM][i])
                        end
                    end
                    ps_data.w = calc_w0(ps_data)
                    ps_data.prim = get_prim(ps_data.w, ka.kinfo)
                    if i==1
                        ps_data.neighbor.state[2] = BALANCE_FLAG
                        ps_data.neighbor.data[2] = [Odata]
                    end
                end
            end
            insert!(datas, index - 1 + i, ps_data)
            dp[] = P4estPsData(pointer_from_objref(ps_data))
        end
    elseif Odata.bound_enc==0 # Outside search radius
        flag = Odata.neighbor.state[1]==BALANCE_FLAG
        for i = 1:(2^DIM)
            pw_in_quad = PointerWrapper(in_quads[i])
            dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
            if flag
                ps_data = Odata.neighbor.data[1][i]
            else
                ps_data = ps_copy(Odata)
                ps_data.ds .*= 0.5
                vs_data = ps_data.vs_data
                @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
                for k in axes(vs_data.df, 2)
                    for j in axes(vs_data.df, 1)
                        @views vs_data.df[j, k]+=min(
                            abs(
                                vs_data.df[j, k]/(
                                    0.5*dot(ps_data.ds, abs.(vs_data.sdf[j, k, :]))+EPS
                                ),
                            ),
                            1.0,
                        )*dot(vs_data.sdf[j, k, :], 0.5 * ps_data.ds .* RMT[DIM][i])
                    end
                end
                ps_data.w = calc_w0(ps_data)
                ps_data.prim = get_prim(ps_data.w, ka.kinfo)
                if i==1
                    ps_data.neighbor.state[2] = BALANCE_FLAG
                    ps_data.neighbor.data[2] = [Odata]
                end
            end
            insert!(datas, index - 1 + i, ps_data)
            dp[] = P4estPsData(pointer_from_objref(ps_data))
        end
    else # Inside search radius
        flag = Odata.neighbor.state[1]==BALANCE_FLAG
        ib = ka.kinfo.config.IB[Odata.bound_enc]
        for i = 1:(2^DIM)
            pw_in_quad = PointerWrapper(in_quads[i])
            dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
            if flag
                ps_data = Odata.neighbor.data[1][i]
            else
                ds = Odata.ds*0.5
                midpoint = @. Odata.midpoint+0.5 * ds * RMT[DIM][i]
                if solid_flag(ib, midpoint) # Inside solid
                    if ghost_cell_flag(ib, midpoint, ds)
                        ps_data = ps_copy(Odata)
                        ps_data.bound_enc = -Odata.bound_enc
                        ps_data.ds .*= 0.5
                        vs_data = ps_data.vs_data
                        @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
                        if i==1
                            ps_data.neighbor.state[2] = BALANCE_FLAG
                            ps_data.neighbor.data[2] = [Odata]
                        end
                    else
                        ps_data = InsideSolidData{DIM,NDF}(-Odata.bound_enc, midpoint, ds)
                    end
                else
                    ps_data = ps_copy(Odata)
                    ps_data.bound_enc = Odata.bound_enc
                    ps_data.ds .*= 0.5
                    vs_data = ps_data.vs_data
                    @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
                    for k in axes(vs_data.df, 2)
                        for j in axes(vs_data.df, 1)
                            @views vs_data.df[j, k]+=min(
                                abs(
                                    vs_data.df[j, k]/(
                                        0.5*dot(ps_data.ds, abs.(vs_data.sdf[j, k, :]))+EPS
                                    ),
                                ),
                                1.0,
                            )*dot(vs_data.sdf[j, k, :], 0.5 * ps_data.ds .* RMT[DIM][i])
                        end
                    end
                    ps_data.w = calc_w0(ps_data)
                    ps_data.prim = get_prim(ps_data.w, ka.kinfo)
                    if i==1
                        ps_data.neighbor.state[2] = BALANCE_FLAG
                        ps_data.neighbor.data[2] = [Odata]
                    end
                end
            end
            insert!(datas, index - 1 + i, ps_data)
            dp[] = P4estPsData(pointer_from_objref(ps_data))
        end
    end
    return nothing
end

function solid_static_replace!(
    ::ChildNum,
    out_quad,
    in_quads,
    which_tree,
    ka::KA{DIM,NDF},
) where {DIM,NDF} # coarsen replace, average or interpolate? Currently interpolation strategy is adopted. If my memory serves me right, problems came out with average most likely due to the iterative balance process.
    trees = ka.kdata.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_in_quad = PointerWrapper(in_quads[1])
    dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
    Odatas = Vector{AbstractPsData{DIM,NDF}}(undef, 2^DIM)
    for i = 1:(2^DIM)
        pw_out_quad = PointerWrapper(out_quad[i])
        Odatas[i] = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4estPsData, pw_out_quad.p.user_data[]).ps_data),
        )
    end
    ds = Odatas[1].ds*2.0
    midpoint = Odatas[1].midpoint-0.5*Odatas[1].ds .* RMT[DIM][1]
    if all(x->isa(x, InsideSolidData), Odatas)
        ps_data = InsideSolidData{DIM,NDF}(Odatas[1].bound_enc, midpoint, ds)
    else
        if Odatas[1].neighbor.state[2]==BALANCE_FLAG
            ps_data = first(Odatas[1].neighbor.data[2])
        else
            _, index = findmin([x.vs_data.vs_num for x in Odatas])
            ps_data = ps_merge(Odatas, index, ka.kinfo)
            vs_data = ps_data.vs_data
            for i in eachindex(Odatas)
                vs_data_n = Odatas[i].vs_data
                vs_merge!(vs_data.sdf, vs_data_n.sdf, vs_data.level, vs_data_n.level, ka)
                vs_merge!(vs_data.df, vs_data_n.df, vs_data.level, vs_data_n.level, ka)
            end
            vs_data.df./=2^DIM;
            vs_data.sdf./=2^DIM
            ps_data.neighbor.state[1] = BALANCE_FLAG
            ps_data.neighbor.data[1] = Odatas
        end
    end
    index = findfirst(x -> x === Odatas[1], datas)
    deleteat!(datas, index:(index+2^DIM-1))
    insert!(datas, index, ps_data)
    dp[] = P4estPsData(pointer_from_objref(ps_data))
    return nothing
end
function initialize_solid_vs_data(ib::AbstractBoundary, kinfo::KInfo)
    prim = get_bc(ib.bc)
    return initialize_vs_data(prim, kinfo)
end

"""
$(SIGNATURES)
Refine grids at original position of the immersed boundaries. The main purpose of the functionality is to obtain
high-resolution field from a fast-obtained coarsener field.
"""
function solid_static_refine!(ka::KA{3})
    kinfo = ka.kinfo
    p4est = kinfo.forest.p4est
    trees = ka.kdata.field.trees.data
    for tree in trees
        for ps_data in tree
            isa(ps_data, InsideSolidData) && continue
            if ps_data.bound_enc<0
                filter_solid_cell!(ps_data, ka)
            end
        end
    end
    c_refine_flag_fn = @cfunction(
        $solid_static_refine_flag,
        Cint,
        (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})
    )
    c_replace_fn = @cfunction(
        $solid_static_replace!,
        Cvoid,
        (
            Ptr{p8est_t},
            p4est_topidx_t,
            Cint,
            Ptr{Ptr{p8est_quadrant_t}},
            Cint,
            Ptr{Ptr{p8est_quadrant_t}},
        )
    )
    GC.@preserve c_refine_flag_fn c_replace_fn p8est_refine_ext(
        p4est,
        1,
        kinfo.config.solver.AMR_PS_MAXLEVEL,
        c_refine_flag_fn,
        C_NULL,
        c_replace_fn,
    )

    c_coarsen_flag_fn = @cfunction(
        $solid_static_coarsen_flag,
        Cint,
        (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}})
    )
    GC.@preserve c_coarsen_flag_fn c_replace_fn p8est_coarsen_ext(
        p4est,
        1,
        0,
        c_coarsen_flag_fn,
        C_NULL,
        c_replace_fn,
    )

    GC.@preserve c_replace_fn p8est_balance_ext(
        p4est,
        P8EST_CONNECT_FULL,
        C_NULL,
        c_replace_fn,
    )
    for tree in trees
        for ps_data in tree
            isa(ps_data, InsideSolidData) && continue
            if ps_data.bound_enc>0
                ps_data.bound_enc = 0
            end
        end
    end
    return nothing
end
