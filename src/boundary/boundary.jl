get_bc(bc::AbstractVector) = bc
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
function colinear_reject!(IB_nodes::Vector{AbstractIBNodes})
    DIM = length(first(IB_nodes).midpoint)
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
function sort_IB_cells!(circle::AbstractCircle,config::Configure,solid_cells::SolidCells,::Vector,IB_cells::IBCells)
    for i in eachindex(solid_cells.ps_datas)
        ps_data = solid_cells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        if isa(config.IB_sort,DistanceIBSort)
            distances = [norm(x.midpoint .-aux_point) for x in IB_cells.IB_nodes[i]]
        end
        index = sortperm(distances)
        IB_cells.IB_nodes[i] = IB_cells.IB_nodes[i][index]
        if isa(config.IB_interp,BilinearIBInterpolate)
            colinear_reject!(IB_cells.IB_nodes[i])
        end
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
        IB_nodes_temp = [AbstractIBNodes[] for _ in solid_cells[i].ps_datas]
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