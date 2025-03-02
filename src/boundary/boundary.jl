get_bc(bc::AbstractVector) = bc
function calc_IB_ρw(aux_point::AbstractVector,bound::Circle,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_2D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Sphere,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_3D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end

function calc_IB_ρw_2D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_2D2F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function calc_IB_ρw_3D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_3D1F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function IB_prim(circle::AbstractCircle,aux_point::AbstractVector,ρw::Real)
    IB_prim(circle.bc,aux_point,ρw)
end
function IB_prim(bc::AbstractVector,::AbstractVector,ρw::Real)
    prim = copy(bc)
    prim[1] = ρw
    return prim
end
function calc_solid_cell_slope!(svdata::AbstractVsData{DIM,NDF},fvdata::VS_Data{DIM,NDF},smid::Vector{Float64},fmid::Vector{Float64},direction::Integer) where{DIM,NDF}
    j = 1
    flag = 0.0
    level = svdata.level
    sdf = @view(svdata.sdf[:,:,direction])
    df = svdata.df
    level_n = fvdata.level
    df_n = fvdata.df
    dx = fmid[direction]-smid[direction]
    @inbounds for i in 1:svdata.vs_num
        if level[i] == level_n[j]
            @views @. sdf[i, :] = df_n[j,:]-df[i,:]
            j += 1
        elseif level[i] < level_n[j]
            @views sdf[i, :] .= -df[i,:]
            while flag != 1.0
                @views @. sdf[i, :] += df_n[j, :]/ 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @views @. sdf[i, :] = df_n[j,:]-df[i,:]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
    sdf./=dx
end
function calc_w0(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,::Global_Data{2,2})
    @views micro_to_macro_2D2F(midpoint[:,1],midpoint[:,2],df[:,1],df[:,2],weight)
end
function calc_qf(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,prim::AbstractVector,::Global_Data{2,2})
    @views heat_flux_2D2F(midpoint[:,1],midpoint[:,2],df[:,1],df[:,2],prim,weight)
end
function calc_pressure(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,::Global_Data{2})
    @views pressure_2D(midpoint[:,1],midpoint[:,2],df[:,1],weight)
end
function save_boundary_result!(ib::Circle,ps_data,solid_neighbor::SolidNeighbor{DIM,NDF,ID},boundary_results,amr::AMR{DIM,NDF}) where{DIM,NDF,ID}
    global_data = amr.global_data
    vs_data = ps_data.vs_data
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    Θ = heaviside.(vn)
    ib_df = vs_data.df
    aux_df = zeros(vs_data.vs_num,NDF)
    dir = get_dir(ID)
    if solid_neighbor.average_num==0
        vs_interpolate!(ib_df,vs_data.level,ps_data.midpoint[dir],s_vs_data.df,
            s_vs_data.level,solid_cell.midpoint[dir],aux_df,aux_point[dir],amr)
        ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
        aux_prim = IB_prim(ib,aux_point,ρw)
        for i in 1:vs_data.vs_num
            if Θ[i]==1.
                aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),aux_prim,amr.global_data)
            end
        end
    else
        for i in 1:vs_data.vs_num
            if Θ[i]==0.
                aux_df[i,:] .= ib_df[i,:]
            end
        end
        ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
        aux_prim = IB_prim(ib,aux_point,ρw)
        for i in 1:vs_data.vs_num
            if Θ[i]==1.
                aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),aux_prim,amr.global_data)
            end
        end
    end
    aux_w = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    aux_prim = get_prim(aux_w,global_data)
    aux_qf = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,aux_prim,global_data)
    aux_p = calc_pressure(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    push!(boundary_results[ps_data.bound_enc].midpoints,aux_point)
    push!(boundary_results[ps_data.bound_enc].ps_solutions,Boundary_PS_Solution(aux_prim,aux_qf,aux_p))
end
function save_boundary_result!(ib::AbstractBoundary,ps_data::PS_Data{DIM,NDF},boundary_results::Vector{Boundary_Solution},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        save_boundary_result!(ib,ps_data,ps_data.neighbor.data[i][1],boundary_results,amr)
    end
end
function solid_cell_index_encoder!(solid_cell_index::Vector{Int},now_index::Int)
    id = findfirst(x->x==0,solid_cell_index)
    isnothing(id) && (@error `A larger SOLID_CELL_ID_NUM is needed!`)
    solid_cell_index[id]=now_index
end
function solid_cell_index_decoder(solid_cell_index::Vector{Int})
    ids = findall(x->x!=0,solid_cell_index)
    return solid_cell_index[ids]
end
function solid_cell_index2ranks(indices::Vector,quadids::Vector,gfq::Vector) # give the 1-based ranks containing solid_cells in indices
    qids = quadids[indices]
    lrank = MPI.Comm_rank(MPI.COMM_WORLD)+1
    ranks = Int[];lids = Int[]
    for i in eachindex(qids)
        rank = findfirst(x->x>qids[i],gfq)-1
        !(rank==lrank||in(rank,ranks))&&push!(ranks,rank)
        rank==lrank&&push!(lids,indices[i])
    end
    return ranks,lids
end
function solid_flag(boundary::AbstractCircle,midpoint::AbstractVector) # Does midpoint locate at solid?
    return xor(norm(midpoint.-boundary.center)>boundary.radius,boundary.solid)
end
function solid_flag(::Domain{Maxwellian},::AbstractVector) # Does midpoint locate at solid?
    return false
end
function solid_cell_flag(boundary::AbstractCircle,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data,inside::Bool) # Ghost nodes, those are inside solid domain and immediately adjacent the boundary.
    (boundary_flag(boundary,midpoint,ds,global_data) && xor(boundary.solid,!inside)) && return true
    return false
end
function solid_cell_flag(boundaries::Vector{AbstractBoundary},midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data)
    for boundary in boundaries
        solid_cell_flag(boundary,midpoint,ds,global_data,solid_flag(boundary,midpoint))&& return true
    end
    return false
end
function InsideSolid_flag(boundary::AbstractBoundary,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data)
    inside = solid_flag(boundary,midpoint)
    (solid_flag(boundary,midpoint)&&!solid_cell_flag(boundary,midpoint,ds,global_data,inside))&& return true # In solid region and not the ghost cell( a.k.a. solid cell)
    return false
end
function InsideSolid_flag(boundaries::Vector{AbstractBoundary},midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data)
    for boundary in boundaries
        InsideSolid_flag(boundary,midpoint,ds,global_data) && return true
    end
    return false
end

function calc_intersect_point(boundaries::Vector{AbstractBoundary},solid_midpoints::Vector)
    aux_points = Vector{Vector{Vector{Float64}}}(undef,length(boundaries))
    for i in eachindex(solid_midpoints)
        aux_points[i] = calc_intersect_point(boundaries[i],solid_midpoints[i])
    end
    return aux_points
end
function calc_intersect_point(boundary::AbstractBoundary,solid_midpoints::Vector{Vector{T}})where{T<:Real}
    aux_points = Vector{Vector{Float64}}(undef,length(solid_midpoints))
    for i in eachindex(solid_midpoints)
        aux_points[i] = calc_intersect_point(boundary,solid_midpoints[i])
    end
    return aux_points
end
function calc_intersect_point(boundary::AbstractCircle,midpoint::AbstractVector{Float64})
    boundary.radius/norm(midpoint-boundary.center)*(midpoint-boundary.center)+boundary.center
end
function calc_normal(midpoint::AbstractVector,circle::Circle)
    aux_point = calc_intersect_point(circle,midpoint)
    return (aux_point-circle.center)/circle.radius # outer normal direction
end
function pre_broadcast_boundary_points(boundary_points::Vector)
    Nb = length(boundary_points)
    numbers = Vector{Int}(undef,Nb)
    for i in eachindex(boundary_points)
        numbers[i] = length(boundary_points[i])
    end
    Numbers = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    for i in eachindex(Numbers)
        if i-1 == MPI.Comm_rank(MPI.COMM_WORLD)
            Numbers[i] = numbers
        else
            Numbers[i] = Vector{Int}(undef,Nb)
        end
    end
    for i in eachindex(Numbers)
        MPI.Bcast!(Numbers[i],i-1,MPI.COMM_WORLD)
    end
    return Numbers
end
function broadcast_boundary_midpoints!(boundary_points::Vector{Vector{Vector{Float64}}},image_points::Vector{Vector{Vector{Float64}}},::Global_Data{DIM})where{DIM}
    Numbers = pre_broadcast_boundary_points(boundary_points)
    rbuffer = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Float64}}(undef,length(boundary_points)) # boundaries{points}
    for i in eachindex(boundary_points)
        buffer = Vector{Float64}(undef,2*DIM*length(boundary_points[i]))
        for j in eachindex(boundary_points[i])
            buffer[2*DIM*(j-1)+1:2*DIM*j-DIM] .= boundary_points[i][j]
            buffer[2*DIM*(j-1)+DIM+1:2*DIM*j] .= image_points[i][j]
        end
        sbuffer[i] = buffer
    end
    for i in eachindex(boundary_points)
        rbuffer[i] = Vector{Vector{Float64}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
        for j in eachindex(Numbers)
            if j-1==MPI.Comm_rank(MPI.COMM_WORLD) 
                rbuffer[i][j] = sbuffer[i]
            else
                rbuffer[i][j] = Vector{Float64}(undef,2*DIM*Numbers[j][i])
            end
        end
    end
    for i in eachindex(boundary_points)
        for j in eachindex(Numbers)
            MPI.Bcast!(rbuffer[i][j],j-1,MPI.COMM_WORLD)
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
    aux_points_global = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points))
    image_points_global = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points))
    for i in eachindex(boundary_points)
        aux_points_global[i] = Vector{Float64}[]
        image_points_global[i] = Vector{Float64}[]
        for j in eachindex(Numbers)
            if j-1 == MPI.Comm_rank(MPI.COMM_WORLD)
                append!(aux_points_global[i],boundary_points[i])
                append!(image_points_global[i],image_points[i])
            else
                for k in 1:Int(length(rbuffer[i][j])/DIM/2)
                    push!(aux_points_global[i],rbuffer[i][j][2*DIM*(k-1)+1:2*DIM*(k-1)+DIM]) 
                    push!(image_points_global[i],rbuffer[i][j][2*DIM*(k-1)+DIM+1:2*DIM*k])
                end
            end
        end
    end
    return aux_points_global,image_points_global
end

function pre_broadcast_quadid(solid_cells::Vector{T}) where{T<:SolidCells}
    Nb = length(solid_cells)
    numbers = Vector{Int}(undef,Nb)
    for i in eachindex(solid_cells)
        numbers[i] = length(solid_cells[i].ps_datas)
    end
    Numbers = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    for i in eachindex(Numbers)
        if i-1 == MPI.Comm_rank(MPI.COMM_WORLD)
            Numbers[i] = numbers
        else
            Numbers[i] = Vector{Int}(undef,Nb)
        end
    end
    for i in eachindex(Numbers)
        MPI.Bcast!(Numbers[i],i-1,MPI.COMM_WORLD)
    end
    return Numbers
end
function upwind_weight(points,aux_point,prim)
    # sum(@views [dot(aux_point-x,prim[2:end-1])/norm((x-aux_point)) for x in points[2:end]])
    average_point = sum(points)/length(points)
    sum([sum((p-average_point).^2) for p in points])*dot(aux_point-average_point,prim[2:end-1])
end
function upwind_template(templates,ib_nodes,aux_point,prim)
    weight = -Inf;template = Int[]
    for temp in templates
        w = upwind_weight([x.midpoint for x in ib_nodes[temp]],aux_point,prim)
        if w>weight
            template = temp;weight = w
        end
    end
    # if aux_point[1]>0&&aux_point[2]>0
    #     @show aux_point [x.midpoint for x in ib_nodes[template]] prim
    # end
    # isnothing(template)&&throw(`No upwind template found. Larger search coefficient is requested.`)
    return template,weight
end
function broadcast_quadid!(solid_cells::Vector{T})where{T<:SolidCells}
    Numbers = pre_broadcast_quadid(solid_cells)
    rbuffer = Vector{Vector{Vector{Cint}}}(undef,length(solid_cells)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Cint}}(undef,length(solid_cells)) # boundaries{points}
    for i in eachindex(solid_cells)
        buffer = Vector{Cint}(undef,length(solid_cells[i].ps_datas))
        for j in eachindex(solid_cells[i].ps_datas)
            buffer[j] = solid_cells[i].quadids[j]
        end
        sbuffer[i] = buffer
    end
    for i in eachindex(solid_cells)
        rbuffer[i] = Vector{Vector{Cint}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
        for j in eachindex(Numbers)
            if j-1==MPI.Comm_rank(MPI.COMM_WORLD) 
                rbuffer[i][j] = sbuffer[i]
            else
                rbuffer[i][j] = Vector{Cint}(undef,Numbers[j][i])
            end
        end
    end
    for i in eachindex(solid_cells)
        for j in eachindex(Numbers)
            MPI.Bcast!(rbuffer[i][j],j-1,MPI.COMM_WORLD)
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
    for i in eachindex(solid_cells)
        for j in eachindex(Numbers)
            j-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
            append!(solid_cells[i].quadids,rbuffer[i][j])
        end
        sort!(solid_cells[i].quadids)
    end
    return Numbers
end

function IB_flag(boundary::AbstractCircle,aux_point::AbstractVector,midpoint::AbstractVector,::AbstractVector)
    r = boundary.search_radius
    norm(midpoint-aux_point)<r
end
function IB_flag(boundaries::Vector{AbstractBoundary},aux_points::Vector{Vector{Vector{Float64}}},midpoint::AbstractVector,ds::AbstractVector)
    for i in eachindex(boundaries)
        solid_flag(boundaries[i],midpoint) && return false
        for j in eachindex(aux_points[i])
            IB_flag(boundaries[i],aux_points[i][j],midpoint,ds) && return true
        end
    end
    return false
end
function calc_image_point(solid_midpoint::AbstractVector{Float64},aux_point::AbstractVector{Float64})
    return 2*aux_point-solid_midpoint
end
function calc_image_point(solid_midpoints::AbstractVector,aux_points::AbstractVector)
    image_points = Vector{Vector{Vector{Float64}}}(undef,length(solid_midpoints))
    for i in eachindex(image_points)
        image_points[i] = Vector{Vector{Float64}}(undef,length(solid_midpoints[i]))
        for j in eachindex(image_points[i])
            image_points[i][j] = calc_image_point(solid_midpoints[i][j],aux_points[i][j])
        end
    end
    return image_points
end

function update_IB_sdata!(boundary::Boundary)
    sdata = boundary.IB_buffer.sdata
    IB_ranks = boundary.IB_ranks
    for i in eachindex(IB_ranks)
        i-1==MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(IB_ranks[i])&&continue
        offset = 0
        for j in eachindex(IB_ranks[i])
            sdata[i].prims[j,:] .= IB_ranks[i][j].prim
            vs_num = IB_ranks[i][j].vs_data.vs_num
            sdata[i].df[offset+1:offset+vs_num,:] .= IB_ranks[i][j].vs_data.df
            offset+=vs_num
        end
    end
end
function IB_data_exchange!(amr::AMR)
    boundary = amr.field.boundary
	update_IB_sdata!(boundary)
    IB_ranks = boundary.IB_ranks
    sdata = boundary.IB_buffer.sdata
    rdata = boundary.IB_buffer.rdata
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sdata)
        i-1 == rank&&continue
        isempty(IB_ranks[i])&&continue
        sreq = MPI.Isend(
            sdata[i].prims,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + rank,
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].df,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + rank,
        )
        push!(reqs, sreq)
    end
    for i in eachindex(rdata)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        !isassigned(r_vs_nums,i) && continue
        rreq = MPI.Irecv!(
                rdata[i].prims,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                rdata[i].df,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
end
function IB_vs_nums_exchange!(boundary::Boundary)
    IB_ranks = boundary.IB_ranks
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    sbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{solid_cells{vs_data{}}}
    for i in eachindex(sbuffer)
        i-1 == rank&&continue
        sbuffer[i] = [x.vs_data.vs_num for x in IB_ranks[i]]
    end
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sbuffer)
        i-1 == rank&&continue
        isempty(sbuffer[i])&&continue
        sreq = MPI.Isend(
            sbuffer[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    for i in eachindex(r_vs_nums)
        i-1 == rank&&continue
        !isassigned(r_vs_nums,i) && continue
        rreq = MPI.Irecv!(
                r_vs_nums[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        MPI.Waitall!(reqs)
    end
end
function IB_structure_exchange!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    boundary = amr.field.boundary
    IB_ranks = boundary.IB_ranks
    sdata = boundary.IB_buffer.sdata
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    reqs = Vector{MPI.Request}(undef, 0)
    vnums = zeros(Int,length(IB_ranks))
    for i in eachindex(IB_ranks)
        i-1==rank&&continue
        isempty(IB_ranks[i])&&continue
        for j in eachindex(IB_ranks[i])
            vnums[i]+=IB_ranks[i][j].vs_data.vs_num
        end
        vs_nums = vnums[i]
        vs_midpoint = Matrix{Cdouble}(undef,vs_nums,DIM)
        df = Matrix{Cdouble}(undef,vs_nums,NDF)
        level = Vector{Int8}(undef,vs_nums)
        sdata[i].level = level
        sdata[i].vs_midpoint = vs_midpoint
        sdata[i].df = df
        offset = 0
        for j in eachindex(IB_ranks[i])
            vs_num = IB_ranks[i][j].vs_data.vs_num
            df[offset+1:offset+vs_num,:] .= IB_ranks[i][j].vs_data.df
            level[offset+1:offset+vs_num] .= IB_ranks[i][j].vs_data.level
            vs_midpoint[offset+1:offset+vs_num,:] .= IB_ranks[i][j].vs_data.midpoint
            offset+=vs_num
        end
        sreq = MPI.Isend(
            sdata[i].level,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].vs_midpoint,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].df,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    rdata = boundary.IB_buffer.rdata
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    for i in eachindex(rdata)
        i-1 == rank&&continue
        !isassigned(r_vs_nums,i) && continue
        vs_nums = sum(r_vs_nums[i])
        rdata[i].df = df = Matrix{Cdouble}(undef,vs_nums,NDF)
        rdata[i].vs_midpoint = vs_midpoint = Matrix{Float64}(undef,vs_nums,DIM)
        rdata[i].level = level = Vector{Int8}(undef,vs_nums)
        rreq = MPI.Irecv!(
                level,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                vs_midpoint,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                df,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
end
function IB_wrap_update!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    boundary = amr.field.boundary
    quadrature = amr.global_data.config.quadrature
    vs_trees_num = reduce(*, amr.global_data.config.vs_trees_num)
    vs_space = 1.0
    for i = 1:DIM
        vs_space *= quadrature[2*i] - quadrature[2*i-1]
    end
    max_weight = vs_space / vs_trees_num
    rdata = boundary.IB_buffer.rdata
    ghost_nodes = boundary.IB_buffer.ghost_nodes
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(rdata)
        i-1==rank&&continue
        !isassigned(r_vs_nums,i)&&continue
        vs_offset = 0
        for j in eachindex(r_vs_nums[i])
            vs_num = r_vs_nums[i][j]
            level = @views rdata[i].level[vs_offset+1:vs_offset+vs_num]
            weight = max_weight./2.0.^(DIM*(level))
            ghost_nodes[i][j].vs_data = @views GhostIBVSData{DIM,NDF}(vs_num,level,weight,
            rdata[i].vs_midpoint[vs_offset+1:vs_offset+vs_num,:],rdata[i].df[vs_offset+1:vs_offset+vs_num,:])
            vs_offset+=vs_num
        end
    end
end

function IB_structure_update!(amr::AMR)
    boundary = amr.field.boundary
    IB_vs_nums_exchange!(boundary)
    IB_structure_exchange!(amr)
    IB_wrap_update!(amr)
    reinit_solid_cells!(boundary)
end
function update_quadid!_kernel(ip,data,dp)
    quadid = global_quadid(ip)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    isa(ps_data,InsideSolidData)&&return nothing
    ps_data.quadid = quadid
    return nothing
end
function update_quadid!(info,data)
    AMR_volume_iterate(info, data, P4est_PS_Data, update_quadid!_kernel)
end
function IB_quadid_exchange!(boundary::Boundary)
    rank = MPI.Comm_rank(MPI.COMM_WORLD) 
    solid_cells = boundary.solid_cells
    Numbers = boundary.solid_numbers
    rbuffer = Vector{Vector{Vector{Cint}}}(undef,length(solid_cells)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Cint}}(undef,length(solid_cells)) # boundaries{points}
    for i in eachindex(solid_cells)
        buffer = [ps_data.quadid for ps_data in solid_cells[i].ps_datas]
        solid_cells[i].quadids = sbuffer[i] = buffer
    end
    for i in eachindex(solid_cells)
        rbuffer[i] = Vector{Vector{Cint}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
        for j in eachindex(Numbers)
            if j-1==rank
                rbuffer[i][j] = sbuffer[i]
            else
                rbuffer[i][j] = Vector{Cint}(undef,Numbers[j][i])
            end
        end
    end
    for i in eachindex(solid_cells)
        for j in eachindex(Numbers)
            isempty(rbuffer[i][j]) && continue
            MPI.Bcast!(rbuffer[i][j],j-1,MPI.COMM_WORLD)
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
    for i in eachindex(solid_cells)
        for j in eachindex(Numbers)
            j-1 == rank&&continue
            append!(solid_cells[i].quadids,rbuffer[i][j])
        end
        sort!(solid_cells[i].quadids)
    end
end
function IB_quadid_update!(ps4est::P_pxest_t,amr::AMR) # Before partition, global quadids of solid_cells and IB_nodes need to be updated
    boundary = amr.field.boundary
    AMR_4est_volume_iterate(ps4est, C_NULL, update_quadid!)
    IB_quadid_exchange!(boundary)
end
function IB_solid_reassign!(ps4est::P_pxest_t,amr::AMR{DIM,NDF}) where{DIM,NDF} # collect IB_nodes for each solid cells, and update local solid cells
    pp = PointerWrapper(ps4est)
    gfq = Base.unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    solid_cells = amr.field.boundary.solid_cells
    Numbers = amr.field.boundary.solid_numbers
    for i in eachindex(solid_cells)
        solid_cell = solid_cells[i]
        solid_cell.ps_datas = PS_Data{DIM,NDF}[]
        for j in 1:MPI.Comm_size(MPI.COMM_WORLD)
            lid = findfirst(x->x>=gfq[j]&&x<gfq[j+1],solid_cell.quadids)
            rid = findfirst(x->x>=gfq[j+1],solid_cell.quadids)
            if isnothing(lid)
                Numbers[j][i]=0
            elseif isnothing(rid)
                Numbers[j][i] = length(solid_cell.quadids)-lid+1
                # Numbers[j+1:end][i] .= 0
                for k in j+1:length(Numbers)
                    Numbers[k][i] = 0
                end
                break
            else
                Numbers[j][i] = rid-lid
            end
        end
    end
    IB_nodes = [PS_Data{DIM,NDF}[] for _ in 1:length(solid_cells)]
    trees = amr.field.trees.data
    for i in eachindex(trees)
        for j in eachindex(trees[i])
            ps_data = trees[i][j]
            isa(ps_data,InsideSolidData)&&continue
            if ps_data.bound_enc>0
                push!(IB_nodes[ps_data.bound_enc],ps_data)
            elseif ps_data.bound_enc<0
                push!(solid_cells[-ps_data.bound_enc].ps_datas,ps_data)
            end
        end
    end
    return IB_nodes
end
function IB_partition!(ps4est::P_pxest_t,amr::AMR{DIM,NDF}) where{DIM,NDF} # After partition, IB need to be reinitialized
    IB_nodes = IB_solid_reassign!(ps4est,amr)
    boundary = amr.field.boundary
    solid_cells = boundary.solid_cells
    IB_ranks,IB_ranks_table,IB_cells = IB_Numbers_ranks(ps4est,IB_nodes,solid_cells)
    boundary.IB_buffer = init_IBCells(ps4est,solid_cells,IB_ranks,IB_ranks_table,IB_cells,amr.global_data)
    boundary.IB_ranks = IB_ranks;boundary.IB_cells = IB_cells
    sort_IB_cells!(amr.global_data,boundary)
    reinit_solid_cells!(boundary)
end
function colinear_test(p1,p2,p3)
    n2 = sum(p1.^2-p2.^2)
    isapprox((p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]))/n2,0.;atol=1e-3)
end
function colinear_reject!(IB_nodes::Vector{AbstractIBNodes{DIM,NDF}}) where{DIM,NDF}
    if DIM==2
        index = zeros(Int,4)
        index[1],index[2] = 1,2
        for i in 3:length(IB_nodes)
            if !colinear_test(IB_nodes[1].midpoint,IB_nodes[2].midpoint,IB_nodes[i].midpoint)
                index[3] = i
                for j in i+1:length(IB_nodes)
                    for k in 1:3,l in (k+1):3
                        if colinear_test(IB_nodes[index[k]].midpoint,IB_nodes[index[l]].midpoint,IB_nodes[j].midpoint)
                            index[4] = 0
                            break
                        end
                        index[4] = j
                    end
                    index[4]!=0&&break
                end
            end
            index[4]!=0 && break
        end
        index[4]==0&& (@error `All nodes combinations are colinear!`)
        IB_nodes[3],IB_nodes[index[3]] = IB_nodes[index[3]],IB_nodes[3]
        IB_nodes[4],IB_nodes[index[4]] = IB_nodes[index[4]],IB_nodes[4]
    elseif DIM==3
    end
end
function singular_test(p1,p2,p3,p4,::BilinearIBInterpolate)
    A = bilinear_coeffi_2D(p1,p2,p3,p4)
    rank(A)<4&&return true
    # @show det(A)
    return false
end
function available_templates!(IB_cell::IBCells{DIM},interp_t::AbstractIBInterpolateType) where{DIM}
    IB_nodes = IB_cell.IB_nodes
    for solid_id in eachindex(IB_cell.IB_nodes)
        IB_nodes = IB_cell.IB_nodes[solid_id]
        if DIM==2
            templates = Vector{Int}[]
            for i in 2:length(IB_nodes), j in (i+1):length(IB_nodes), k in (j+1):length(IB_nodes)
                index = collect(1:4)
                index[2:4] .= i,j,k
                # for l in 1:4, m in l+1:4, n in m+1:4
                    # if singular_test(IB_nodes[index[l]].midpoint,IB_nodes[index[m]].midpoint,IB_nodes[index[n]].midpoint,interp_t)
                    # if singular_test(IB_nodes[index[1]].midpoint,IB_nodes[index[2]].midpoint,IB_nodes[index[3]].midpoint,IB_nodes[index[4]].midpoint,interp_t)
                    #     flag = false
                    # end
                # end
                if !singular_test(IB_nodes[index[1]].midpoint,IB_nodes[index[2]].midpoint,IB_nodes[index[3]].midpoint,IB_nodes[index[4]].midpoint,interp_t)
                    push!(templates,index) 
                end
            end
            length(templates)<1&& throw(`Search coefficient is requested to be larger.`)
            IB_cell.templates[solid_id] = templates
        elseif DIM==3
        end
    end
end
function sort_IB_cells!(circle::AbstractCircle,config::Configure,solid_cells::SolidCells,::Vector,IB_cells::IBCells)
    for i in eachindex(solid_cells.ps_datas)
        ps_data = solid_cells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        if isa(config.IB_sort,DistanceIBSort)
            distances = [norm(x.midpoint .-aux_point) for x in IB_cells.IB_nodes[i]]
        end
        index = sortperm(distances)
        IB_cells.IB_nodes[i] = IB_cells.IB_nodes[i][index]
        # if isa(config.IB_interp,BilinearIBInterpolate)
        #     colinear_reject!(IB_cells.IB_nodes[i])
        #     # available_templates!(IB_cells)
        # end
        available_templates!(IB_cells,config.IB_interp)
    end
end
function sort_IB_cells!(global_data::Global_Data,boundary::Boundary)
    IBs = global_data.config.IB
    solid_cells = boundary.solid_cells
    aux_points = boundary.aux_points
    IB_cells = boundary.IB_cells
    for i in eachindex(IBs)
        sort_IB_cells!(IBs[i],global_data.config,solid_cells[i],aux_points[i],IB_cells[i])
    end
end

function init_solid_cells!(boundary::Boundary{DIM,NDF},global_data::Global_Data{DIM,NDF})where{DIM,NDF}
    IBs = global_data.config.IB
    solid_cells = boundary.solid_cells
    IB_cells = boundary.IB_cells
    @inbounds for i in eachindex(solid_cells)
        for j in eachindex(solid_cells[i].ps_datas)
            ps_data = solid_cells[i].ps_datas[j]
            vs_data = ps_data.vs_data
            IB_node = first(IB_cells[i].IB_nodes[j])
            IB_vs = IB_node.vs_data
            vs_data.vs_num = IB_vs.vs_num
            vs_data.level = copy(IB_vs.level)
            vs_data.weight = copy(IB_vs.weight)
			vs_data.midpoint = copy(IB_vs.midpoint)
            ps_data.prim = get_bc(IBs[i].bc)
            vs_data.df = discrete_maxwell(ps_data,global_data)
			vs_data.sdf = zeros(IB_vs.vs_num,NDF,DIM) 
            vs_data.flux = zeros(IB_vs.vs_num,NDF)
        end
    end
end
function reinit_solid_cells!(boundary::Boundary{DIM,NDF})where{DIM,NDF}
    solid_cells = boundary.solid_cells
    IB_cells = boundary.IB_cells
    @inbounds for i in eachindex(solid_cells)
        for j in eachindex(solid_cells[i].ps_datas)
            ps_data = solid_cells[i].ps_datas[j]
            vs_data = ps_data.vs_data
            IB_node = first(IB_cells[i].IB_nodes[j])
            IB_vs = IB_node.vs_data
            df_temp = zeros(IB_vs.vs_num,NDF)
            vs_projection!(IB_vs,vs_data,df_temp)
            vs_data.df = df_temp
            vs_data.vs_num = IB_vs.vs_num
            vs_data.level = copy(IB_vs.level)
            vs_data.weight = copy(IB_vs.weight)
			vs_data.midpoint = copy(IB_vs.midpoint)
			vs_data.sdf = zeros(vs_data.vs_num,NDF,DIM) 
        end
    end
end


function init_IB!(ps4est::P_pxest_t,trees::PS_Trees{DIM,NDF},global_data::Global_Data{DIM,NDF},solid_cells::Vector{T},aux_points::Vector) where{DIM,NDF,T<:SolidCells}
    solid_Numbers = broadcast_quadid!(solid_cells)
    IB_nodes = [PS_Data{DIM,NDF}[] for _ in 1:length(solid_cells)]
    search_IB!(IB_nodes,aux_points,trees,global_data)
    IB_ranks,IB_ranks_table,IB_cells = IB_Numbers_ranks(ps4est,IB_nodes,solid_cells)
    IB_buffer = init_IBCells(ps4est,solid_cells,IB_ranks,IB_ranks_table,IB_cells,global_data)
    return solid_Numbers,IB_cells,IB_buffer,IB_ranks
end
function search_IB!(IB_nodes::Vector{Vector{PS_Data{DIM,NDF}}},aux_points::Vector,trees::PS_Trees{DIM,NDF},global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    ps_datas = trees.data
    boundaries = global_data.config.IB
    for i in eachindex(ps_datas)
        for j in eachindex(ps_datas[i])
            ps_data = ps_datas[i][j]
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0) && continue
            for k in eachindex(boundaries)
                for l in eachindex(aux_points[k])
                    if IB_flag(boundaries[k],aux_points[k][l],ps_data.midpoint,ps_data.ds)
                        ps_data.bound_enc = k
                        first(ps_data.solid_cell_index)==0&&push!(IB_nodes[k],ps_data)
                        solid_cell_index_encoder!(ps_data.solid_cell_index,l)
                    end
                end
            end
        end
    end
end
function IB_Numbers_ranks(ps4est::P_pxest_t,IB_nodes::Vector,solid_cells::Vector{SolidCells{DIM,NDF}}) where{DIM,NDF}
    pp = PointerWrapper(ps4est)
    gfq = unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    IB_ranks_table = [[PS_Data{DIM,NDF}[] for _ in 1:length(solid_cells)] for _ in 1:length(gfq)-1]
    IB_cells = Vector{IBCells}(undef,length(solid_cells))
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(solid_cells)
        local_offset = findfirst(x->x>=gfq[rank+1]&&x<gfq[rank+2],solid_cells[i].quadids)
        IB_nodes_temp = [AbstractIBNodes{DIM,NDF}[] for _ in solid_cells[i].ps_datas]
        for IB_node in IB_nodes[i]
            ids = solid_cell_index_decoder(IB_node.solid_cell_index)
            ranks,lids = solid_cell_index2ranks(ids,solid_cells[i].quadids,gfq)
            for l in lids
                push!(IB_nodes_temp[l-local_offset+1],IB_node)
            end
            for r in ranks
                push!(IB_ranks_table[r][i],IB_node)
            end
        end
        IB_cells[i] = IBCells(IB_nodes_temp)
    end
    IB_ranks = [[x for u in v for x in u] for v in IB_ranks_table]
    return IB_ranks,IB_ranks_table,IB_cells
end
function init_IBCells(ps4est::P_pxest_t,solid_cells::Vector{T},IB_ranks::Vector,IB_ranks_table::Vector,IB_cells::Vector{IBCells},global_data::Global_Data) where{T<:SolidCells}
    rNumbers,r_vs_nums,IB_buffer = IB_nodes_communicate(solid_cells,IB_ranks,IB_ranks_table) # rNumbers: Vector{Vector{Int}}
    ghost_nodes = IB_nodes_data_wrap!(ps4est,rNumbers,r_vs_nums,IB_buffer.rdata,solid_cells,IB_cells,global_data)
    IB_buffer.r_vs_nums = r_vs_nums;IB_buffer.ghost_nodes = ghost_nodes
    return IB_buffer
end
function IB_nodes_communicate(solid_cells::Vector{T},IB_ranks::Vector,IB_ranks_table::Vector) where{T<:SolidCells}
    rNumbers,s_vs_nums,r_vs_nums = IB_nodes_pre_communicate(solid_cells,IB_ranks_table) # Communicate Numbers and vs_nums for buffer allocation and topological info quadids.
    IB_buffer = IB_nodes_data_communicate(IB_ranks,s_vs_nums,r_vs_nums)
    return rNumbers,r_vs_nums,IB_buffer
end
function IB_nodes_pre_communicate(solid_cells::Vector{T},IB_ranks_table::Vector) where{T<:SolidCells}
    rNumbers = IB_nodes_Numbers_communicate(IB_ranks_table,solid_cells)
    s_vs_nums,r_vs_nums = IB_nodes_vs_nums_communicate(rNumbers,IB_ranks_table)
    return rNumbers,s_vs_nums,r_vs_nums
end
function IB_nodes_Numbers_communicate(IB_ranks_table::Vector,solid_cells::Vector{T})where{T<:SolidCells}
    rbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{boundaries{}}
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(rbuffer)
        i-1==rank&&continue
        rbuffer[i] = Vector{Int}(undef,length(solid_cells))
    end
    sbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    for i in eachindex(rbuffer)
        i-1 == rank && continue
        sbuffer[i] = [length(x) for x in IB_ranks_table[i]]
    end
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sbuffer)
        i-1 == rank&&continue
        sreq = MPI.Isend(
            sbuffer[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    for i in eachindex(rbuffer)
        i-1 == rank&&continue
        rreq = MPI.Irecv!(
                rbuffer[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
    return rbuffer
end
function IB_nodes_vs_nums_communicate(rNumbers::Vector,IB_ranks_table::Vector)
    sbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{boundaries{IB_nodes{}}}
    for i in eachindex(sbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        sbuffer[i] = [x.vs_data.vs_num for v in IB_ranks_table[i] for x in v]
    end
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(sbuffer[i])&&continue
        sreq = MPI.Isend(
            sbuffer[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    r_vs_nums = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # ranks{IB_nodes{}}
    for i in eachindex(r_vs_nums)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        (num = sum(rNumbers[i]))==0 && continue
        r_vs_nums[i] = Vector{Int}(undef,num)
        rreq = MPI.Irecv!(
                r_vs_nums[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
    return sbuffer,r_vs_nums
end
function IB_nodes_data_communicate(IB_ranks::Vector{Vector{PS_Data{DIM,NDF}}},s_vs_nums::Vector,r_vs_nums::Vector) where{DIM,NDF}
    IB_buffer = collect_IB_nodes_data(IB_ranks,s_vs_nums)
    sdata = IB_buffer.sdata
    IB_buffer.rdata = rdata = Vector{IBTransferData}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sdata)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(IB_ranks[i])&&continue
        sreq = MPI.Isend(
            sdata[i].solid_cell_indices,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].midpoints,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].prims,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].level,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].vs_midpoint,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].df,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    for i in eachindex(rdata)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        !isassigned(r_vs_nums,i) && continue
        ps_num = length(r_vs_nums[i])
        vs_nums = sum(r_vs_nums[i])
        solid_cell_indices = Matrix{Int}(undef,ps_num,SOLID_CELL_ID_NUM)
        midpoints = Matrix{Cdouble}(undef,ps_num,DIM)
        prims = Matrix{Cdouble}(undef,ps_num,DIM+2)
        df = Matrix{Cdouble}(undef,vs_nums,NDF)
        vs_midpoint = Matrix{Cdouble}(undef,vs_nums,DIM)
        level = Vector{Int8}(undef,vs_nums)
        rdata[i] = IBTransferData(solid_cell_indices,midpoints,prims,level,vs_midpoint,df)
        rreq = MPI.Irecv!(
                solid_cell_indices,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                midpoints,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                prims,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                level,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                vs_midpoint,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
                df,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
    return IB_buffer
end
function collect_IB_nodes_data(IB_ranks::Vector{Vector{PS_Data{DIM,NDF}}},s_vs_nums::Vector) where{DIM,NDF}
    sdata = Vector{IBTransferData}(undef,length(IB_ranks)) # ranks{solid_cells{}}, data distributes as [df(sum(vs_nums),NDF),level(sum(vs_nums))]
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    IB_buffer = IBBuffer();
    for i in eachindex(IB_ranks)
        i-1==rank&&continue
        isempty(IB_ranks[i])&&continue
        vs_nums = sum(s_vs_nums[i])
        ps_num = length(IB_ranks[i])
        solid_cell_indices = Matrix{Int}(undef,ps_num,SOLID_CELL_ID_NUM)
        midpoints = Matrix{Cdouble}(undef,ps_num,DIM)
        prims = Matrix{Cdouble}(undef,ps_num,DIM+2)
        df = Matrix{Cdouble}(undef,vs_nums,NDF)
        vs_midpoint = Matrix{Cdouble}(undef,vs_nums,DIM)
        level = Vector{Int8}(undef,vs_nums)
        sdata[i] = IBTransferData(solid_cell_indices,midpoints,prims,level,vs_midpoint,df)
        offset = 0
        for j in eachindex(IB_ranks[i])
            solid_cell_indices[j,:] .= IB_ranks[i][j].solid_cell_index
            midpoints[j,:] .= IB_ranks[i][j].midpoint
            prims[j,:] .= IB_ranks[i][j].prim
            vs_num = IB_ranks[i][j].vs_data.vs_num
            df[offset+1:offset+vs_num,:] .= IB_ranks[i][j].vs_data.df
            vs_midpoint[offset+1:offset+vs_num,:] .= IB_ranks[i][j].vs_data.midpoint
            level[offset+1:offset+vs_num] .= IB_ranks[i][j].vs_data.level
            offset+=vs_num
        end
    end
    IB_buffer.sdata = sdata
    return IB_buffer
end
function IB_nodes_data_wrap!(ps4est::P_pxest_t,rNumbers::Vector{Vector{Int}},r_vs_nums::Vector{Vector{Int}},rdata::Vector{IBTransferData},solid_cells::Vector{SolidCells{DIM,NDF}},IB_cells::Vector{IBCells},global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    quadrature = global_data.config.quadrature
    vs_trees_num = reduce(*, global_data.config.vs_trees_num)
    vs_space = 1.0
    for i = 1:DIM
        vs_space *= quadrature[2*i] - quadrature[2*i-1]
    end
    max_weight = vs_space / vs_trees_num
    pp = PointerWrapper(ps4est)
    gfq = unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    ghost_nodes = Vector{Vector{GhostIBNode{DIM,NDF}}}(undef,length(rdata))
    for i in eachindex(rdata)
        i-1==rank&&continue
        !isassigned(r_vs_nums,i)&&continue
        ps_offset = 1;vs_offset = 0
        ghost_nodes[i] = Vector{GhostIBNode{DIM,NDF}}(undef,sum(rNumbers[i]))
        for j in eachindex(IB_cells)
            local_offset = findfirst(x->x>=gfq[rank+1]&&x<gfq[rank+2],solid_cells[j].quadids)
            sup = findlast(x->x>=gfq[rank+1]&&x<gfq[rank+2],solid_cells[j].quadids)
            for _ in 1:rNumbers[i][j]
                vs_num = r_vs_nums[i][ps_offset]
                level = @views rdata[i].level[vs_offset+1:vs_offset+vs_num]
                weight = max_weight./2.0.^(DIM*level)
                ghost_nodes[i][ps_offset] = ghost_node = @views GhostIBNode{DIM,NDF}(
                    rdata[i].solid_cell_indices[ps_offset,:],rdata[i].midpoints[ps_offset,:],
                    rdata[i].prims[ps_offset,:],GhostIBVSData{DIM,NDF}(vs_num,level,
                    weight,rdata[i].vs_midpoint[vs_offset+1:vs_offset+vs_num,:],
                    rdata[i].df[vs_offset+1:vs_offset+vs_num,:]))
                solid_cell_indices = solid_cell_index_decoder(ghost_node.solid_cell_index)
                ls = findall(x->x>=local_offset&&x<=sup,solid_cell_indices)
                local_ids = solid_cell_indices[ls].-local_offset.+1
                for local_id in local_ids
                    push!(IB_cells[j].IB_nodes[local_id],ghost_node)
                end
                ps_offset+=1;vs_offset+=vs_num
            end
        end 
    end
    return ghost_nodes
end



function vs_extrapolate!(df::AbstractMatrix{Float64},sdf::AbstractMatrix{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},dx::Float64,::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(dft,1)
        if levelt[i] == level[j]
            @. dft[i, :] = @views df[j, :]+sdf[j,:]*dx
            j += 1
        elseif levelt[i] < level[j]
            while flag != 1.0
                @. dft[i, :] += @views (df[j, :]+sdf[j,:]*dx)/ 2^(DIM * (level[j] - levelt[i]))
                flag += 1 / 2^(DIM * (level[j] - levelt[i]))
                j += 1
            end
            flag = 0.0
        else
            @. dft[i, :] = @views df[j,:]+sdf[j,:]*dx
            flag += 1 / 2^(DIM * (levelt[i] - level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function vs_extrapolate!(df::AbstractMatrix{Float64},sdf::AbstractMatrix{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},dx::Float64,weight::AbstractVector{Float64},::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(dft,1)
        if levelt[i] == level[j]
            @. dft[i, :] += @views (df[j, :]+sdf[j,:]*dx)*weight[i]
            j += 1
        elseif levelt[i] < level[j]
            while flag != 1.0
                @. dft[i, :] += @views (df[j, :]+sdf[j,:]*dx)/ 2^(DIM * (level[j] - levelt[i]))*weight[i]
                flag += 1 / 2^(DIM * (level[j] - levelt[i]))
                j += 1
            end
            flag = 0.0
        else
            @. dft[i, :] += @views (df[j,:]+sdf[j,:]*dx)*weight[i]
            flag += 1 / 2^(DIM * (levelt[i] - level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function update_solid_cell!(::Val{1},ps_data::PS_Data{DIM,NDF},fluid_cells::Vector,dirs::Vector{Int},::Vector{Float64},amr::AMR{DIM,NDF}) where{DIM,NDF}
    vs_data = ps_data.vs_data
    f_vs_data = fluid_cells[1].vs_data
    dir = dirs[1]
    fdf = f_vs_data.df;fsdf = @views f_vs_data.sdf[:,:,dir];dx = ps_data.midpoint[dir]-fluid_cells[1].midpoint[dir]
    vs_data.df.=0.
    vs_extrapolate!(fdf,fsdf,f_vs_data.level,vs_data.df,vs_data.level,dx,amr)    
end
function update_solid_cell!(::Val{2},ps_data::PS_Data{2,NDF},fluid_cells::Vector,dirs::Vector{Int},rots::Vector{Float64},amr::AMR{2,NDF}) where{NDF}
    if dirs[1]==dirs[2]
        throw(`The solid region is too thin that only one solid_cell is contained.`)
    end
    vs_data = ps_data.vs_data;vs_data.df.=0.
    # weight = [u[1]^2/sum(u.^2) for u in eachrow(vs_data.midpoint)] # cosθ^2
    weightx = [abs(u[1])/norm(u) for u in eachrow(vs_data.midpoint)] # cosθ
    weighty = [abs(u[2])/norm(u) for u in eachrow(vs_data.midpoint)] # sinθ
    weights = Matrix{Float64}(undef,length(weightx),length(fluid_cells))
    for i in eachindex(fluid_cells)
        dir = dirs[i];rot = rots[i]
        # weights[:,i] .= (dir==1 ? copy(weight) : 1 .-weight)
        weights[:,i] .= (dir==1 ? copy(weightx) : copy(weighty))
        opp_id = findall(x->x*rot>0,@views vs_data.midpoint[:,dir])
        weights[opp_id,i] .= 0.
    end
    weight_i = Vector{Float64}(undef,length(weightx))
    weight_sum = sum(weights,dims=2)
    for i in eachindex(fluid_cells)
        f_vs_data = fluid_cells[i].vs_data
        dir = dirs[i]
        for j in eachindex(weight_i)
            weight_i[j] = weight_sum[j]==0 ? 1/length(fluid_cells) : weights[j,i]/weight_sum[j]
        end
        fdf = f_vs_data.df;fsdf = @views f_vs_data.sdf[:,:,dir];dx = ps_data.midpoint[dir]-fluid_cells[i].midpoint[dir]
        vs_extrapolate!(fdf,fsdf,f_vs_data.level,vs_data.df,vs_data.level,dx,weight_i,amr)
    end
end
function update_solid_cell!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = [ps_data.neighbor.data[i][1] for i in fluid_dirs]
    dirs = get_dir.(fluid_dirs)
    rots = get_rot.(fluid_dirs)
    update_solid_cell!(Val(length(fluid_dirs)),ps_data,fluid_cells,dirs,rots,amr)
end
function update_solid_cell!(amr::AMR)
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc>=0)&&continue
            update_solid_cell!(ps_data,amr)
        end
    end
end
function initialize_solid_neighbor!(amr::AMR)
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            initialize_solid_neighbor!(ps_data,amr)
        end
    end
end
function calc_intersect(f_midpoint,s_midpoint,circle::Circle)
    c = circle.center;r = circle.radius
    if abs(f_midpoint[1]-s_midpoint[1])<EPS
        t = acos((f_midpoint[1]-c[1])/r)
        if f_midpoint[2]>c[2]
            ap = [r*cos(t),r*sin(t)];n=ap./r
        else
            ap = [r*cos(t),-r*sin(t)];n = ap./r
        end
    else
        t = asin((f_midpoint[2]-c[2])/r)
        if f_midpoint[1]>c[1]
            ap = [r*cos(t),r*sin(t)];n = ap./r
        else
            ap = [-r*cos(t),r*sin(t)];n = ap./r
        end
    end
    return ap,n
end
function initialize_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_dirs = findall(x->!isnothing(x[1])&&x[1].bound_enc<0,ps_data.neighbor.data)
    vs_data = ps_data.vs_data
    for i in solid_dirs
        solid_cell = ps_data.neighbor.data[i][1]
        if ps_data.bound_enc!=0
            ps_data.bound_enc!=-solid_cell.bound_enc&&throw(`Different immersed boundaries correspond to the same fluid cell!`)
        end
        ps_data.bound_enc=-solid_cell.bound_enc
        ib = amr.global_data.config.IB[-solid_cell.bound_enc]
        aux_point,normal = calc_intersect(ps_data.midpoint,solid_cell.midpoint,ib)
        dir = get_dir(i)
        av_num = abs(solid_cell.midpoint[dir]-ps_data.midpoint[dir])<10*ps_data.ds[dir]^2 ? 1 : 0
        ic = get_bc(ib.bc)
        svsdata = VS_Data{DIM,NDF}(
            vs_data.vs_num,
            vs_data.level,
            vs_data.weight,
            vs_data.midpoint,
            discrete_maxwell(vs_data.midpoint,ic,amr.global_data),
            zeros(vs_data.vs_num,NDF,DIM),
            # zeros(vs_data.vs_num,NDF,DIM),
            Matrix{Float64}(undef,0,0)
        )
        ps_data.neighbor.data[i][1] = SolidNeighbor{DIM,NDF,i}(
            solid_cell.bound_enc,av_num,aux_point,normal,solid_cell,solid_cell.midpoint,svsdata
        )
        # if any(x->isnan(x),ps_data.neighbor.data[i][1].vs_data.df)
        #     throw(`initial sn.df nan!`)
        # end
    end
    av_id = findall(x->ps_data.neighbor.data[x][1].average_num==1,solid_dirs)
    for id in av_id
        i = solid_dir[id]
        ps_data.neighbor.data[i][1].average_num = length(av_id)
    end
end
function vs_interpolate!(f_df::AbstractMatrix,f_level::AbstractVector{Int8},fx,s_df,s_level,sx,b_df,bx,::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(f_df,1)
        if f_level[i] == s_level[j]
            @. b_df[i, :] = @views f_df[i,:]+(s_df[j, :]-f_df[i,:])/(sx-fx)*(bx-fx)
            j += 1
        elseif f_level[i] < s_level[j]
            @. b_df[i,:] = @views -f_df[i,:]
            while flag != 1.0
                @. b_df[i, :] += @views s_df[j,:]/ 2^(DIM * (s_level[j] - f_level[i]))
                flag += 1 / 2^(DIM * (s_level[j] - f_level[i]))
                j += 1
            end
            @. b_df[i,:]*=(bx-fx)/(sx-fx)
            @. b_df[i,:] += @views f_df[i,:]
            flag = 0.0
        else
            @. b_df[i, :] += @views f_df[i,:]+(s_df[j, :]-f_df[i,:])/(sx-fx)*(bx-fx)
            flag += 1 / 2^(DIM * (f_level[i] - s_level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF,ID},amr::AMR) where{DIM,NDF,ID}
    ib = amr.global_data.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    # dxL = norm(aux_point-ps_data.midpoint)
    dir = get_dir(ID)
    dxL = ps_data.ds[dir]
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    dxR = norm(solid_cell.midpoint-aux_point)
    ib_df = vs_data.df
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    Θ = heaviside.(vn)
    if solid_neighbor.average_num==0
        vs_interpolate!(ib_df,vs_data.level,ps_data.midpoint[dir],s_vs_data.df,
            s_vs_data.level,solid_cell.midpoint[dir],aux_df,aux_point[dir],amr)
        ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
        aux_prim = IB_prim(ib,aux_point,ρw)
        for i in 1:vs_data.vs_num
            if Θ[i]==1.
                aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),aux_prim,amr.global_data)
            end
        end
        @. solid_neighbor.vs_data.df = aux_df+(aux_df-vs_data.df)/dxL*dxR
    else
        for i in 1:vs_data.vs_num
            if Θ[i]==0.
                aux_df[i,:] .= ib_df[i,:]
            end
        end
        ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
        aux_prim = IB_prim(ib,aux_point,ρw)
        for i in 1:vs_data.vs_num
            if Θ[i]==1.
                aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),aux_prim,amr.global_data)
            end
        end
        @. solid_neighbor.vs_data.df = aux_df
    end
    # end
end
function vs_average!(dfs::Vector,levels,df,level,::AMR{DIM,NDF}) where{DIM,NDF}
    num = length(dfs)
    for k in eachindex(dfs)
        j = 1;flag = 0.0
        dfn = dfs[k]
        leveln = levels[k]
        @inbounds for i in axes(df,1)
            if level[i] == leveln[j]
                @. df[i, :] = @views dfn[j,:]/num
                j += 1
            elseif level[i] < leveln[j]
                while flag != 1.0
                    @. df[i, :] += @views dfn[j,:]/ 2^(DIM * (leveln[j] - level[i]))/num
                    flag += 1 / 2^(DIM * (level[j] - leveln[i]))
                    j += 1
                end
                flag = 0.0
            else
                @. df[i, :] += @views dfn[j,:]/num
                flag += 1 / 2^(DIM * (level[i] - leveln[j]))
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
        end
    end
end
function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        update_solid_neighbor!(ps_data,ps_data.neighbor.data[i][1],amr)
    end
    av_ids = findall(x->ps_data.neighbor.data[x][1].average_num!=0,solid_neighbors)
    if !isempty(av_ids)
        # num = length(av_ids)
        ps_data.vs_data.df .= 0.
        dfs = [ps_data.neighbor.data[solid_neighbors[i]][1].vs_data.df for i in av_ids]
        levels = [ps_data.neighbor.data[solid_neighbors[i]][1].vs_data.level for i in av_ids]
        # for i in av_ids
        #     @. ps_data.vs_data.df+=ps_data.neighbor.data[solid_neighbors[i]][1].vs_data.df/num
        # end
        vs_average!(dfs,levels,ps_data.vs_data.df,ps_data.vs_data.level,amr)
    end
end
function update_solid_neighbor!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
            update_solid_neighbor!(ps_data,amr)
        end
    end
end
function update_solid!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    initialize_solid_neighbor!(amr)
end