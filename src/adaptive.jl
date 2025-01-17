# function boundary_flag(amr::AMR{DIM,NDF}, ps_data::PS_Data{DIM,NDF}) where{DIM,NDF}
#     geometry = amr.global_data.config.geometry
#     for i = 1:DIM
#         abs(ps_data.midpoint[i] - geometry[2*i-1]) < ps_data.ds[i] && return true
#         abs(ps_data.midpoint[i] - geometry[2*i]) < ps_data.ds[i] && return true
#     end
#     return false
# end
function coarsen_vs_average!(sdf::AbstractArray,sdf_n::AbstractArray,level::Vector,level_n::Vector,::AMR{DIM}) where{DIM}
    j = 1
    flag = 0.0
    for i in axes(sdf,1)
        if level[i] == level_n[j]
            @. sdf[i,:,:] += @views sdf_n[j, :, :]
            j += 1
        elseif level[i] < level_n[j]
            while flag != 1.0
                @. sdf[i, :, :] += @views sdf_n[j, :, :]/ 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @. sdf[i, :, :] += @views sdf_n[j,:,:]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function coarsen_vs_average!(sdf::AbstractMatrix,sdf_n::AbstractMatrix,level::Vector,level_n::Vector,::AMR{DIM}) where{DIM}
    j = 1
    flag = 0.0
    for i in axes(sdf,1)
        if level[i] == level_n[j]
            @. sdf[i,:] += @views sdf_n[j, :]
            j += 1
        elseif level[i] < level_n[j]
            while flag != 1.0
                @. sdf[i, :] += @views sdf_n[j, :]/ 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @. sdf[i, :] += @views sdf_n[j,:]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function boundary_flag(boundary::Domain,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data) # Domain boundary flag
    index = Int(floor((boundary.id-1)/2))+1
    abs(midpoint[index] - global_data.config.geometry[boundary.id]) < ds[index] && return true
end
function boundary_flag(boundary::Circle,midpoint::AbstractVector,ds::AbstractVector,::Global_Data) # Circle type IB boundary flag
    #dsmin = [global_data.config.geometry[2i]-global_data.config.geometry[2i-1] for i in 1:2]/2^global_data.config.solver.AMR_PS_MAXLEVEL./global_data.config.trees_num
    flag = 0
    for i = 1:4
        #flag += norm(midpoint.+0.5*(ds+dsmin).*RMT[2][i].-boundary.center)>boundary.radius ? 1 : -1 # any corner cross boundary?
        flag += norm(midpoint.+ds.*NMT[2][i].-boundary.center)>boundary.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==4 && return false
    return true
end
function domain_flag(global_data::Global_Data,midpoint::AbstractVector,ds::AbstractVector) # domain flag for physical space dynamic adaptive
    domains = global_data.config.domain
    for domain in domains
        boundary_flag(domain,midpoint,ds,global_data) && return true
    end
    return false
end
function boundary_flag(global_data::Global_Data,midpoint::AbstractVector,ds::AbstractVector)
    domain = global_data.config.domain
    boundary = global_data.config.IB
    for i in eachindex(domain)
        boundary_flag(domain[i],midpoint,ds,global_data) && return true
    end
    return solid_cell_flag(boundary,midpoint,ds,global_data)
end
function boundary_flag(
    fp::PW_pxest_t,
    treeid,
    qp::PW_pxest_quadrant_t,
    global_data::Global_Data{DIM},
) where{DIM}
    domain = global_data.config.domain
    boundary = global_data.config.IB
    ds, midpoint = quad_to_cell(fp, treeid, qp)
    for i in eachindex(domain)
        boundary_flag(domain[i],midpoint,ds,global_data) && return Cint(1)
    end
    for i in eachindex(boundary)
        boundary_flag(boundary[i],midpoint,ds,global_data) && return Cint(1)
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
function ps_merge(Odatas::Vector,index::Int,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    data = Odatas[index]
    vs_data = data.vs_data;vs_num = vs_data.vs_num
    p = PS_Data(DIM,NDF,
        data.ds *2.0,
        data.midpoint - 0.5 * data.ds .* RMT[DIM][index],
        [sum([x.w[i] for x in Odatas]) for i in 1:DIM+2]./=2^DIM,
        [sum([x.sw[i,j] for x in Odatas]) for i in 1:DIM+2,j in 1:DIM]./=2^DIM,
        VS_Data{DIM,NDF}(vs_num,copy(vs_data.level),copy(vs_data.weight),copy(vs_data.midpoint),
            zeros(vs_num,NDF),zeros(vs_num,NDF,DIM),zeros(vs_num,NDF))
    )
    # p.neighbor = data.neighbor
    p.neighbor = Neighbor(DIM,NDF)
    # if !(first(data.neighbor.state)==BALANCE_FLAG)
    #     p.neighbor.state[1] = BALANCE_FLAG
    #     p.neighbor.data[1] = Odatas
    # end
    p.neighbor.state[1] = BALANCE_FLAG
    p.neighbor.data[1] = Odatas
    p.prim = get_prim(p.w,global_data)
    p.qf = Vector{Float64}(undef, DIM)
    p.flux = zeros(DIM + 2)
    return p
end
function ps_refine_flag(
    ps_data::PS_Data{DIM},
    amr::AMR{DIM},
    qp::PW_pxest_quadrant_t,
) where{DIM}
    global_data = amr.global_data
    if ps_data.bound_enc!=0||domain_flag(global_data,ps_data.midpoint,ps_data.ds)
        return Cint(1)
    end
    # agrad = @views abs.((ps_data.sw[end, :]))
    agrad = [maximum(abs.(ps_data.sw[i,:])) for i in 1:DIM+2]
    rgrad = maximum(agrad ./ (global_data.status.gradmax.+EPS))
    if rgrad > 2.0^(DIM+qp.level[] - global_data.config.solver.AMR_PS_MAXLEVEL) * ADAPT_COEFFI_PS
        flag = Cint(1)
    else
        flag = Cint(0)
    end
    flag
end
function ps_refine_flag(
    ::InsideSolidData,
    ::AMR,
    ::PW_pxest_quadrant_t,
)
    Cint(0)
end
function ps_refine_flag(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ps_refine_flag(ps_data, amr, qp)
    end
end
function ps_replace(::Val{1}, out_quad, in_quads, which_tree, amr::AMR{DIM}) where{DIM}# refine replace
    trees = amr.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_out_quad = PointerWrapper(out_quad[1])
    Odata = unsafe_pointer_to_objref(
        pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
    )
    index = findfirst(x -> x === Odata, datas)
    deleteat!(datas, index)
    flag = first(Odata.neighbor.state)==BALANCE_FLAG
    for i = 1:2^DIM
        pw_in_quad = PointerWrapper(in_quads[i])
        dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
        if flag
            ps_data = Odata.neighbor.data[1][i]
        else
            ps_data = ps_copy(Odata)
            ps_data.ds .*= 0.5
            vs_data = ps_data.vs_data
            @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
            for j = 1:DIM
                @. vs_data.df += 0.5 * ps_data.ds[j] * RMT[DIM][i][j] * @view(vs_data.sdf[:, :, j])
                @. ps_data.w += 0.5 * ps_data.ds[j] * RMT[DIM][i][j] * @view(ps_data.sw[:, j])
            end
            # ps_data.neighbor.state[1] = Odata.neighbor.state[1]
            # ps_data.neighbor.data[1] = Odata.neighbor.data[1]
        end
        insert!(datas, index - 1 + i, ps_data)
        dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    end
end
function ps_replace(::ChildNum, out_quad, in_quads, which_tree, amr::AMR{DIM,NDF}) where{DIM,NDF} # coarsen replace, average or interpolate? Currently interpolation strategy is adopted. If my memory serves me right, problems came out with average most likely due to the iterative balance process.
    trees = amr.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_in_quad = PointerWrapper(in_quads[1])
    dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
    Odatas = Vector{AbstractPsData{DIM,NDF}}(undef, 2^DIM)
    # if !any(x->isa(PS_Data,x),Odatas)
    #     ps_data = copy(Odatas[1])
    # else
    for i = 1:2^DIM
        pw_out_quad = PointerWrapper(out_quad[i])
        Odatas[i] = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
        )
    end
    # @views _,index = findmax([maximum(abs.(x.sw[end,:])) for x in Odatas])
    _,index = findmin([x.vs_data.vs_num for x in Odatas])
    ps_data = ps_merge(Odatas,index,amr.global_data)
    vs_data = ps_data.vs_data
    for i in eachindex(Odatas)
        vs_data_n = Odatas[i].vs_data
        coarsen_vs_average!(vs_data.sdf,vs_data_n.sdf,vs_data.level,vs_data_n.level,amr)
        coarsen_vs_average!(vs_data.df,vs_data_n.df,vs_data.level,vs_data_n.level,amr)
    end
    vs_data.df./=2^DIM;vs_data.sdf./=2^DIM
    # end
    index = findfirst(x -> x === Odatas[1], datas)
    deleteat!(datas, index:index+2^DIM-1)
    insert!(datas, index, ps_data)
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
end
# function ps_balance_replace(::Val{1}, out_quad, in_quads, which_tree, amr::AMR{DIM,NDF}) where{DIM,NDF}# refine replace
#     trees = amr.field.trees
#     treeid = Int(which_tree) - trees.offset
#     datas = trees.data[treeid]
#     pw_out_quad = PointerWrapper(out_quad[1])
#     Odata = unsafe_pointer_to_objref(
#         pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
#     )
#     index = findfirst(x -> x === Odata, datas)
#     deleteat!(datas, index)
#     pw_in_quad = PointerWrapper(in_quads[i])
#     dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
#     if !isa(ps_data,InsideSolidData)
#         if first(Odata.neighbor.state)==BALANCE_FLAG
#             for i = 1:2^DIM
#                 ps_data = Odata.neighbor.data[1][i]
#                 dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
#                 insert!(datas, index - 1 + i, ps_data)
#             end
#         else
#             children = Vector{PS_Data{DIM,NDF}}(undef,2^DIM)
#             for i = 1:2^DIM
#                 ps_data = ps_copy(Odata)
#                 ps_data.ds .*= 0.5
#                 vs_data = ps_data.vs_data
#                 @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
#                 for j = 1:DIM
#                     @. vs_data.df += 0.5 * ps_data.ds[j] * RMT[DIM][i][j] * @view(vs_data.sdf[:, :, j])
#                     @. ps_data.w += 0.5 * ps_data.ds[j] * RMT[DIM][i][j] * @view(ps_data.sw[:, j])
#                 end
#                 children[i] = ps_data
#                 ps_data.neighbor.data
#                 insert!(datas, index - 1 + i, ps_data)
#                 dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
#                 Odata
#             end
#         end
#     end
# end
# function ps_balance_replace(::ChildNum, out_quad, in_quads, which_tree, amr::AMR{DIM,NDF}) where{DIM,NDF} # coarsen replace, average or interpolate? Currently interpolation strategy is adopted. If my memory serves me right, problems came out with average most likely due to the iterative balance process.
#     trees = amr.field.trees
#     treeid = Int(which_tree) - trees.offset
#     datas = trees.data[treeid]
#     pw_in_quad = PointerWrapper(in_quads[1])
#     dp = PointerWrapper(P4est_PS_Data, pw_in_quad.p.user_data[])
#     Odatas = Vector{AbstractPsData{DIM,NDF}}(undef, 2^DIM)
#     # if !any(x->isa(PS_Data,x),Odatas)
#     #     ps_data = copy(Odatas[1])
#     # else
#     for i = 1:2^DIM
#         pw_out_quad = PointerWrapper(out_quad[i])
#         Odatas[i] = unsafe_pointer_to_objref(
#             pointer(PointerWrapper(P4est_PS_Data, pw_out_quad.p.user_data[]).ps_data),
#         )
#     end
#     # @views _,index = findmax([maximum(abs.(x.sw[end,:])) for x in Odatas])
#     _,index = findmax([x.vs_data.vs_num for x in Odatas])
#     Odata = Odatas[index]
#     ps_data = ps_copy(Odata)
#     @. ps_data.midpoint -= 0.5 * ps_data.ds * RMT[DIM][index]
#     # vs_data = ps_data.vs_data
#     # for i = 1:DIM
#     #     @. vs_data.df -= 0.5 * ps_data.ds[i] * RMT[DIM][index][i] * @view(vs_data.sdf[:, :, i])
#     #     @. ps_data.w -= 0.5 * ps_data.ds[i] * RMT[DIM][index][i] * @view(ps_data.sw[:, i])
#     # end
#     ps_data.ds .*= 2.0
#     # end
#     index = findfirst(x -> x === Odatas[1], datas)
#     deleteat!(datas, index:index+2^DIM-1)
#     insert!(datas, index, ps_data)
#     dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
# end
function p4est_replace(forest::T1, which_tree, num_out, out_quads::Ptr{T2}, num_in, in_quads) where{T1<:P_pxest_t,T2<:P_pxest_quadrant_t}
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{T2}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{T2}, in_quads, num_in)
        ps_replace(Val(Int(num_out)), out_quads_wrap, in_quads_wrap, which_tree, amr)
        return nothing
    end
end

function ps_refine!(p4est::Ptr{p4est_t},amr::AMR; recursive = 0)
    p4est_refine_ext(
        p4est,
        recursive,
        amr.global_data.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            ps_refine_flag,
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
            ps_refine_flag,
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


function pre_IB_refine_flag(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        global_data,aux_points = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ds,midpoint = quad_to_cell(fp,which_tree,qp)
        IB_flag(global_data.config.IB,aux_points,midpoint,ds) && return Cint(1)
    end
    return Cint(0)
end
function IB_pre_ps_refine!(p4est::Ptr{p4est_t},global_data::Global_Data)
    p4est_refine_ext(
        p4est,
        1,
        global_data.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            pre_IB_refine_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        C_NULL,
    )
end
function pre_ps_coarsen_flag(forest,which_tree,quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 2
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, quadrants, 2^DIM)
        global_data,aux_points = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        for i = 1:2^DIM
            qp = PointerWrapper(quadrants_wrap[i])
            ds,midpoint = quad_to_cell(fp,which_tree,qp)
            (IB_flag(global_data.config.IB,aux_points,midpoint,ds)||boundary_flag(global_data,midpoint,ds))&&return Cint(0)
        end
        return Cint(1)
    end
end
function pre_ps_coarsen!(p4est::Ptr{p4est_t}; recursive = 0)
    p4est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            pre_ps_coarsen_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}})
        ),
        C_NULL,
        C_NULL,
    )
end

function ps_coarsen_flag(ps_datas::Vector{PS_Data}, levels::Vector{Int}, amr::AMR{DIM,NDF}) where{DIM,NDF}
    global_data = amr.global_data
    flag = Cint(1)
    for i = 1:2^DIM
        # (isa(ps_datas[i],InsideSolidData)||ps_datas[i].bound_enc!=0) && continue
        ps_data = ps_datas[i]
        (ps_data.bound_enc!=0||domain_flag(global_data,ps_data.midpoint,ps_data.ds)) && return Cint(0)
        agrad = [maximum(abs.(ps_data.sw[i,:])) for i in 1:DIM+2]
        rgrad = maximum(agrad ./ (global_data.status.gradmax.+EPS))
        if rgrad > 2.0^(DIM+levels[i] - global_data.config.solver.AMR_PS_MAXLEVEL) * ADAPT_COEFFI_PS
            flag = Cint(0)
            return flag
        end
    end
    return flag
end
function ps_coarsen_flag(forest::Ptr{p4est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 2
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, quadrants, 2^DIM)
        levels = Vector{Int}(undef, 2^DIM)
        ps_datas = Vector{PS_Data}(undef, 2^DIM)
        for i = 1:2^DIM
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4est_PS_Data, qp.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data,InsideSolidData)&&return Cint(0)
            levels[i] = qp.level[]
            ps_datas[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
        end
        amr = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ps_coarsen_flag(ps_datas, levels, amr)
    end
end
function ps_coarsen_flag(forest::Ptr{p8est_t}, which_tree, quadrants)
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
        ps_coarsen_flag(ps_datas, levels, amr)
    end
end
function ps_coarsen!(p4est::Ptr{p4est_t}; recursive = 0)
    p4est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            ps_coarsen_flag,
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
            ps_coarsen_flag,
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
function pre_ps_refine_flag(forest::P_pxest_t, which_tree, quadrant)
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
            pre_ps_refine_flag,
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
            pre_ps_refine_flag,
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
    amr.field.faces = Vector{AbstractFace}(undef, 0)
    initialize_faces!(p4est, amr)
end

