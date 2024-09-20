function boundary_flag(amr::AMR{DIM,NDF}, ps_data::PS_Data{DIM,NDF}) where{DIM,NDF}
    geometry = amr.global_data.config.geometry
    for i = 1:DIM
        abs(ps_data.midpoint[i] - geometry[2*i-1]) < ps_data.ds[i] && return true
        abs(ps_data.midpoint[i] - geometry[2*i]) < ps_data.ds[i] && return true
    end
    return false
end
function boundary_flag(
    fp::PW_pxest_t,
    treeid,
    qp::PW_pxest_quadrant_t,
    global_data::Global_Data{DIM},
) where{DIM}
    geometry = global_data.config.geometry
    ds, midpoint = quad_to_cell(fp, treeid, qp)
    for i = 1:DIM
        abs(midpoint[i] - geometry[2*i-1]) < ds[i] && return Cint(1)
        abs(midpoint[i] - geometry[2*i]) < ds[i] && return Cint(1)
    end
    return Cint(0)
end
function ps_copy(data::PS_Data{DIM,NDF}) where{DIM,NDF}
    p = PS_Data(DIM,NDF,
        copy(data.ds),
        copy(data.midpoint),
        copy(data.w),
        copy(data.sw),
        deepcopy(data.vs_data),
    )
    p.neighbor = Neighbor(DIM,NDF)
    p.prim = copy(data.prim)
    p.qf = Vector{Float64}(undef, DIM)
    p.flux = zeros(DIM + 2)
    return p
end
function refine_flag_ps(
    ps_data::PS_Data{DIM},
    amr::AMR{DIM},
    qp::PW_pxest_quadrant_t,
) where{DIM}
    global_data = amr.global_data
    if boundary_flag(amr, ps_data)
        return Cint(1)
    end
    agrad = maximum(abs.(@view(ps_data.sw[end, :])))
    rgrad = agrad / global_data.status.gradmax
    if rgrad > 2.0^(DIM+qp.level[] - global_data.config.solver.AMR_PS_MAXLEVEL) * 0.01
        flag = Cint(1)
    else
        flag = Cint(0)
    end
    flag
end
function p4est_refine_flag(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        refine_flag_ps(ps_data, amr, qp)
    end
end
function replace_ps(::Val{1}, out_quad, in_quads, which_tree, amr::AMR{DIM}) where{DIM}# refine replace
    trees = amr.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_out_quad = PointerWrapper(out_quad[1])
    Odata = unsafe_pointer_to_objref(
        pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
    )
    index = findfirst(x -> x === Odata, datas)
    deleteat!(datas, index)
    for i = 1:2^DIM
        pw_in_quad = PointerWrapper(in_quads[i])
        dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
        ps_data = ps_copy(Odata)
        ps_data.ds .*= 0.5
        vs_data = ps_data.vs_data
        @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
        for j = 1:DIM
            @. vs_data.df += 0.5 * ps_data.ds[j] * RMT[DIM][i][j] * @view(vs_data.sdf[:, :, j])
            @. ps_data.w += 0.5 * ps_data.ds[j] * RMT[DIM][i][j] * @view(ps_data.sw[:, j])
        end
        insert!(datas, index - 1 + i, ps_data)
        dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    end
end
function replace_ps(::ChildNum, out_quad, in_quads, which_tree, amr::AMR{DIM,NDF}) where{DIM,NDF} # coarsen replace, average or interpolate? Currently interpolation strategy is adopted. If my memory serves me right, problems came out with average most likely due to the iterative balance process.
    trees = amr.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_in_quad = PointerWrapper(in_quads[1])
    dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
    Odatas = Vector{PS_Data{DIM,NDF}}(undef, 2^DIM)
    for i = 1:2^DIM
        pw_out_quad = PointerWrapper(out_quad[i])
        Odatas[i] = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
        )
    end
    index = 1
    Odata = Odatas[index]
    ps_data = ps_copy(Odata)
    @. ps_data.midpoint -= 0.5 * ps_data.ds * RMT[DIM][index]
    vs_data = ps_data.vs_data
    for i = 1:DIM
        @. vs_data.df -= 0.5 * ps_data.ds[i] * RMT[DIM][index][i] * @view(vs_data.sdf[:, :, i])
        @. ps_data.w -= 0.5 * ps_data.ds[i] * RMT[DIM][index][i] * @view(ps_data.sw[:, i])
    end
    ps_data.ds .*= 2.0
    index = findfirst(x -> x === Odatas[1], datas)
    deleteat!(datas, index:index+2^DIM-1)
    insert!(datas, index, ps_data)
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
end
function p4est_replace(forest::T1, which_tree, num_out, out_quads::Ptr{T2}, num_in, in_quads) where{T1<:P_pxest_t,T2<:P_pxest_quadrant_t}
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{T2}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{T2}, in_quads, num_in)
        replace_ps(Val(Int(num_out)), out_quads_wrap, in_quads_wrap, which_tree, amr)
        return nothing
    end
end


function ps_refine!(p4est::Ptr{p4est_t},amr::AMR; recursive = 0)
    p4est_refine_ext(
        p4est,
        recursive,
        amr.global_data.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            p4est_refine_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function ps_refine!(p4est::Ptr{p8est_t},amr::AMR; recursive = 0)
    p8est_refine_ext(
        p4est,
        recursive,
        amr.global_data.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            p4est_refine_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
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
    )
end


function coarsen_flag_ps(ps_datas::Vector{PS_Data}, levels::Vector{Int}, amr::AMR{DIM,NDF}) where{DIM,NDF}
    global_data = amr.global_data
    flag = Cint(1)
    for i = 1:2^DIM
        if boundary_flag(amr, ps_datas[i])
            return Cint(0)
        end
        agrad = maximum(abs.(@view(ps_datas[i].sw[end, :])))
        rgrad = agrad / global_data.status.gradmax
        if rgrad > 2.0^(DIM+levels[i] - global_data.config.solver.AMR_PS_MAXLEVEL) * 0.01
            flag = Cint(0)
            return flag
        end
    end
    return flag
end
function p4est_coarsen_flag(forest::Ptr{p4est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 2
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, quadrants, 2^DIM)
        levels = Vector{Int}(undef, 2^DIM)
        ps_datas = Vector{PS_Data}(undef, 2^DIM)
        for i = 1:2^DIM
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
            levels[i] = qp.level[]
            ps_datas[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
        end
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        coarsen_flag_ps(ps_datas, levels, amr)
    end
end
function p4est_coarsen_flag(forest::Ptr{p8est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 3
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p8est_quadrant_t}}, quadrants, 2^DIM)
        levels = Vector{Int}(undef, 2^DIM)
        ps_datas = Vector{PS_Data}(undef, 2^DIM)
        for i = 1:2^DIM
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
            levels[i] = qp.level[]
            ps_datas[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
        end
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        coarsen_flag_ps(ps_datas, levels, amr)
    end
end
function ps_coarsen!(p4est::Ptr{p4est_t}; recursive = 0)
    p4est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            p4est_coarsen_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function ps_coarsen!(p4est::Ptr{p8est_t}; recursive = 0)
    p8est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            p4est_coarsen_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
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
    )
end
function ps_balance!(p4est::Ptr{p4est_t})
    p4est_balance_ext(
        p4est,
        P4EST_CONNECT_FACE,
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function ps_balance!(p4est::Ptr{p8est_t})
    p8est_balance_ext(
        p4est,
        P8EST_CONNECT_FACE,
        C_NULL,
        @cfunction(
            p4est_replace,
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
    )
end
function pre_p4est_refine_flag(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        global_data = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        boundary_flag(fp, which_tree, qp, global_data)
    end
end
function pre_ps_refine!(p4est::Ptr{p4est_t},global_data::Global_Data)
    p4est_refine_ext(
        p4est,
        1,
        global_data.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            pre_p4est_refine_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        C_NULL,
    )
end
function pre_ps_refine!(p4est::Ptr{p8est_t},global_data::Global_Data)
    p8est_refine_ext(
        p4est,
        1,
        global_data.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            pre_p4est_refine_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})
        ),
        C_NULL,
        C_NULL,
    )
end

function pre_ps_balance!(p4est::Ptr{p4est_t})
    p4est_balance_ext(p4est, P4EST_CONNECT_FACE, C_NULL, C_NULL)
end
function pre_ps_balance!(p4est::Ptr{p8est_t})
    p8est_balance_ext(p4est, P8EST_CONNECT_FACE, C_NULL, C_NULL)
end

function update_faces!(p4est::P_pxest_t, amr::AMR)
    amr.field.faces = Vector{Face}(undef, 0)
    initialize_faces!(p4est, amr)
end
