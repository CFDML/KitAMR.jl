function boundary_flag(DVM_data::DVM_Data, ps_data::AbstractPsData)
    geometry = DVM_data.global_data.geometry
    for i = 1:DIM
        abs(ps_data.midpoint[i] - geometry[2*i-1]) < ps_data.ds[i] && return true
        abs(ps_data.midpoint[i] - geometry[2*i]) < ps_data.ds[i] && return true
    end
    return false
end
function boundary_flag(
    fp::PointerWrapper{p4est_t},
    treeid,
    qp::PointerWrapper{p4est_quadrant_t},
    global_data::Global_Data,
)
    geometry = global_data.geometry
    ds, midpoint = get_midpoint_ds(fp, treeid, qp)
    for i = 1:DIM
        abs(midpoint[i] - geometry[2*i-1]) < ds[i] && return Cint(1)
        abs(midpoint[i] - geometry[2*i]) < ds[i] && return Cint(1)
    end
    return Cint(0)
end
function PS_copy(data::PS_Data)
    p = PS_Data(
        copy(data.ds),
        copy(data.midpoint),
        copy(data.w),
        copy(data.sw),
        deepcopy(data.vs_data),
    )
    p.neighbor = Neighbor(Union{AbstractPsData,Nothing})
    p.prim = copy(data.prim)
    p.qf = Vector{Float64}(undef, DIM)
    p.flux = zeros(DIM + 2)
    return p
end
function refine_flag_PS(
    ps_data::PS_Data,
    DVM_data::DVM_Data,
    qp::PointerWrapper{p4est_quadrant_t},
)
    global_data = DVM_data.global_data
    if boundary_flag(DVM_data, ps_data)
        return Cint(1)
    end
    # agrad = maximum(abs.(ps_data.sw))
    agrad = maximum(abs.(@view(ps_data.sw[4, :])))
    rgrad = agrad / global_data.gradmax
    # if ((rgrad/0.2 > qp.level[]+1)&&(agrad*maximum(ps_data.ds)^2>0.00001))
    # if rgrad/(0.3/DVM_PS_MAXLEVEL) > qp.level[]+1
    if rgrad > 4.0^(qp.level[] - DVM_PS_MAXLEVEL) * 0.01
        flag = Cint(1)
    else
        flag = Cint(0)
    end
    flag
end
function p4est_refine_flag(forest::Ptr{p4est_t}, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        DVM_data = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        refine_flag_PS(ps_data, DVM_data, qp)
    end
end
function replace_PS(::Val{1}, out_quad, in_quads, which_tree, DVM_data::DVM_Data) # refine replace
    treeid = Int(which_tree) - DVM_data.trees.offset
    # @show treeid
    datas = DVM_data.trees.data[treeid]
    pw_out_quad = PointerWrapper(out_quad[1])
    Odata = unsafe_pointer_to_objref(
        pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
    )
    index = findfirst(x -> x === Odata, datas)
    deleteat!(datas, index)
    for i = 1:(2^DIM)
        pw_in_quad = PointerWrapper(in_quads[i])
        dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
        ps_data = PS_copy(Odata)
        ps_data.ds .*= 0.5
        vs_data = ps_data.vs_data
        @. ps_data.midpoint += 0.5 * ps_data.ds * rmt[i]
        for j = 1:DIM
            @. vs_data.df += 0.5 * ps_data.ds[j] * rmt[i][j] * @view(vs_data.sdf[:, :, j])
            @. ps_data.w += 0.5 * ps_data.ds[j] * rmt[i][j] * @view(ps_data.sw[:, j])
        end
        insert!(datas, index - 1 + i, ps_data)
        dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    end
end
function replace_PS(::Val{2^DIM}, out_quad, in_quads, which_tree, DVM_data::DVM_Data) # coarsen replace
    treeid = Int(which_tree) - DVM_data.trees.offset
    datas = DVM_data.trees.data[treeid]
    pw_in_quad = PointerWrapper(in_quads[1])
    dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
    Odatas = Vector{PS_Data}(undef, 2^DIM)
    for i = 1:(2^DIM)
        pw_out_quad = PointerWrapper(out_quad[i])
        Odatas[i] = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
        )
    end
    # _,index = findmax(x->x.vs_data.vs_num,Odatas)
    index = 1
    Odata = Odatas[index]
    ps_data = PS_copy(Odata)
    @. ps_data.midpoint -= 0.5 * ps_data.ds * rmt[index]
    vs_data = ps_data.vs_data
    for i = 1:DIM
        @. vs_data.df -= 0.5 * ps_data.ds[i] * rmt[index][i] * @view(vs_data.sdf[:, :, i])
        @. ps_data.w -= 0.5 * ps_data.ds[i] * rmt[index][i] * @view(ps_data.sw[:, i])
    end
    ps_data.ds .*= 2.0
    index = findfirst(x -> x === Odatas[1], datas)
    deleteat!(datas, index:(index+2^DIM-1))
    insert!(datas, index, ps_data)
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
end
function p4est_replace(forest, which_tree, num_out, out_quads, num_in, in_quads)
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        DVM_data = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, in_quads, num_in)
        replace_PS(Val(Int(num_out)), out_quads_wrap, in_quads_wrap, which_tree, DVM_data)
        return nothing
    end
end

function PS_refine!(p4est::Ptr{p4est_t}; recursive = 0)
    p4est_refine_ext(
        p4est,
        recursive,
        DVM_PS_MAXLEVEL,
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


function coarsen_flag_PS(ps_datas::Vector{PS_Data}, levels::Vector{Int}, DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    flag = Cint(1)
    for i = 1:(2^DIM)
        if boundary_flag(DVM_data, ps_datas[i])
            return Cint(0)
        end
        # agrad = maximum(abs.(ps_datas[i].sw))
        agrad = maximum(abs.(@view(ps_datas[i].sw[4, :])))
        rgrad = agrad / global_data.gradmax
        # if (rgrad/0.2 > levels[i])&&(agrad*maximum(ps_datas[i].ds)^2>0.00001/4.)
        # if rgrad/(0.3/DVM_PS_MAXLEVEL) > levels[i]
        if rgrad > 4.0^(levels[i] - DVM_PS_MAXLEVEL) * 0.01
            flag = Cint(0)
            return flag
        end
    end
    return flag
end
function p4est_coarsen_flag(forest::Ptr{p4est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, quadrants, 2^DIM)
        levels = Vector{Int}(undef, 2^DIM)
        ps_datas = Vector{PS_Data}(undef, 2^DIM)
        for i = 1:(2^DIM)
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
            levels[i] = qp.level[]
            ps_datas[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
        end
        DVM_data = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        coarsen_flag_PS(ps_datas, levels, DVM_data)
    end
end
function PS_coarsen!(p4est::Ptr{p4est_t}; recursive = 0)
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
function PS_balance!(p4est::Ptr{p4est_t})
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
function pre_p4est_refine_flag(forest::Ptr{p4est_t}, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        global_data = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        boundary_flag(fp, which_tree, qp, global_data)
    end
end
function pre_PS_refine!(p4est::Ptr{p4est_t})
    p4est_refine_ext(
        p4est,
        1,
        DVM_PS_MAXLEVEL,
        @cfunction(
            pre_p4est_refine_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        C_NULL,
    )
end
function pre_PS_balance!(p4est::Ptr{p4est_t})
    p4est_balance_ext(p4est, P4EST_CONNECT_FACE, C_NULL, C_NULL)
end
function update_faces!(p4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    DVM_data.faces = Vector{Face}(undef, 0)
    initialize_faces!(p4est, DVM_data)
end
