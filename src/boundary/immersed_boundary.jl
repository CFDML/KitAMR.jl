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
function broadcast_boundary_points!(boundary_points::Vector{Vector{Vector{Float64}}},::Global_Data{DIM})where{DIM}
    Numbers = pre_broadcast_boundary_points(boundary_points)
    rbuffer = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Float64}}(undef,length(boundary_points)) # boundaries{points}
    for i in eachindex(boundary_points)
        buffer = Vector{Float64}(undef,DIM*length(boundary_points[i]))
        for j in eachindex(boundary_points[i])
            buffer[DIM*(j-1)+1:DIM*j] .= boundary_points[i][j]
            # buffer[2*DIM*(j-1)+DIM+1:2*DIM*j] .= normals[i][j]
        end
        sbuffer[i] = buffer
    end
    for i in eachindex(boundary_points)
        rbuffer[i] = Vector{Vector{Float64}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
        for j in eachindex(Numbers)
            if j-1==MPI.Comm_rank(MPI.COMM_WORLD) 
                rbuffer[i][j] = sbuffer[i]
            else
                rbuffer[i][j] = Vector{Float64}(undef,DIM*Numbers[j][i])
            end
        end
    end
    for i in eachindex(boundary_points)
        for j in eachindex(Numbers)
            MPI.Bcast!(rbuffer[i][j],j-1,MPI.COMM_WORLD)
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
    boundary_points_global = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points))
    # normals_global = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points))
    for i in eachindex(boundary_points)
        boundary_points_global[i] = Vector{Float64}[]
        # normals_global[i] = Vector{Float64}[]
        for j in eachindex(Numbers)
            if j-1 == MPI.Comm_rank(MPI.COMM_WORLD)
                append!(boundary_points_global[i],boundary_points[i])
                # append!(normals_global[i],normals[i])
            else
                for k in 1:Int(length(rbuffer[i][j])/DIM)
                    push!(boundary_points_global[i],rbuffer[i][j][DIM*(k-1)+1:DIM*k]) 
                    # push!(normals_global[i],rbuffer[i][j][2*DIM*(k-1)+DIM+1:2*DIM*k])
                end
            end
        end
    end
    return boundary_points_global
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

function IB_node_flag(boundary::AbstractBoundary,aux_point::AbstractVector,midpoint::AbstractVector,::AbstractVector)
    r = boundary.search_radius
    norm(midpoint-aux_point)<r
end
# function IB_flag(boundaries::Vector{AbstractBoundary},boundary_points::Vector{Vector{Vector{Float64}}},midpoint::AbstractVector,ds::AbstractVector)
#     for i in eachindex(boundaries)
#         solid_flag(boundaries[i],midpoint) && return false
#         for j in eachindex(boundary_points[i])
#             IB_flag(boundaries[i],boundary_points[i][j],midpoint,ds) && return true
#         end
#     end
#     return false
# end

function update_IB_sdata!(boundary::ImmersedBoundary)
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
    boundary = amr.field.immersed_boundary
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
    MPI.Waitall(reqs)
end
function IB_vs_nums_exchange!(boundary::ImmersedBoundary)
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
        MPI.Waitall(reqs)
    end
end
function IB_structure_exchange!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    boundary = amr.field.immersed_boundary
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
    MPI.Waitall(reqs)
end
function IB_wrap_update!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    boundary = amr.field.immersed_boundary
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
function IB_quadid_exchange!(boundary::ImmersedBoundary)
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
            if ps_data.IB_nodes_communicate>0
                push!(IB_nodes[ps_data.ib_enc],ps_data)
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
    boundary.IB_buffer = init_IB_cells(ps4est,solid_cells,IB_ranks,IB_ranks_table,IB_cells,amr.global_data)
    boundary.IB_ranks = IB_ranks;boundary.IB_cells = IB_cells
    sort_IB_cells!(amr.global_data,boundary)
    reinit_solid_cells!(boundary)
end

function initialize_solidcells(amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_cells = Vector{SolidCells{DIM,NDF}}(undef,length(amr.global_data.config.IB))
    for i in eachindex(solid_cells)
        solid_cells[i] = SolidCells{DIM,NDF}(PS_Data{DIM,NDF}[],Int[])
    end
    trees = amr.field.trees
    for tree in trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc>=0)&&continue
            push!(solid_cells[-ps_data.bound_enc].ps_datas,ps_data)
            push!(solid_cells[-ps_data.bound_enc].quadids,ps_data.quadid)
        end
    end
    return solid_cells
end
function initialize_IB_nodes(amr::AMR{DIM,NDF}) where{DIM,NDF}
    trees = amr.field.trees
    solid_midpoints_global = amr.field.immersed_boundary.solid_midpoints_global
    boundaries = amr.global_data.config.IB
    IB_nodes = [PS_Data{DIM,NDF}[] for _ in 1:length(boundaries)]
    for tree in trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0) && continue
            for i in eachindex(solid_midpoints_global)
                for j in eachindex(solid_midpoints_global[i])
                    solid_midpoint = solid_midpoints_global[i][j]
                    if IB_node_flag(boundaries[i],solid_midpoint,ps_data.midpoint,ps_data.ds)
                        ps_data.ib_enc = i
                        first(ps_data.solid_cell_index)==0&&push!(IB_nodes[i],ps_data)
                        solid_cell_index_encoder!(ps_data.solid_cell_index,j)
                    end
                end
            end
        end
    end
    return IB_nodes
end
function initialize_IB_ranks(solid_cells,IB_nodes_origin,amr::AMR{DIM,NDF}) where{DIM,NDF}
    pp = PointerWrapper(amr.global_data.forest.p4est)
    gfq = unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    IB_ranks_table = [[PS_Data{DIM,NDF}[] for _ in 1:length(solid_cells)] for _ in 1:length(gfq)-1]
    # IB_cells = Vector{IBCells}(undef,length(solid_cells))
    IB_nodes = Vector{Vector{Vector{AbstractIBNode}}}(undef,length(solid_cells))
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(solid_cells)
        local_offset = findfirst(x->x>=gfq[rank+1]&&x<gfq[rank+2],solid_cells[i].quadids)
        IB_nodes_temp = [AbstractIBNode[] for _ in solid_cells[i].ps_datas]
        for IB_node in IB_nodes_origin[i]
            ids = solid_cell_index_decoder(IB_node.solid_cell_index)
            ranks,lids = solid_cell_index2ranks(ids,solid_cells[i].quadids,gfq)
            for l in lids
                push!(IB_nodes_temp[l-local_offset+1],IB_node)
            end
            for r in ranks
                push!(IB_ranks_table[r][i],IB_node)
            end
        end
        IB_nodes[i] = IB_nodes_temp
    end
    IB_ranks = [[x for u in v for x in u] for v in IB_ranks_table]
    return IB_ranks,IB_ranks_table,IB_nodes
end
function initialize_IB_buffer(solid_cells::Vector{T},IB_ranks::Vector,IB_ranks_table::Vector,IB_nodes::Vector{Vector{Vector{AbstractIBNode}}},amr) where{T<:SolidCells}
    ps4est = amr.global_data.forest.p4est
    rNumbers,r_vs_nums,IB_buffer = IB_nodes_communicate(solid_cells,IB_ranks,IB_ranks_table) # rNumbers: Vector{Vector{Int}}
    ghost_nodes = IB_nodes_data_wrap!(ps4est,rNumbers,r_vs_nums,IB_buffer.rdata,solid_cells,IB_nodes,amr.global_data)
    IB_buffer.r_vs_nums = r_vs_nums;IB_buffer.ghost_nodes = ghost_nodes
    return IB_buffer
end
function initialize_downwind_transport(solid_cells,amr)
    solid_numbers = broadcast_quadid!(solid_cells)
    IB_nodes_origin = initialize_IB_nodes(amr)
    IB_ranks,IB_ranks_table,IB_nodes = initialize_IB_ranks(solid_cells,IB_nodes_origin,amr)
    IB_buffer = initialize_IB_buffer(solid_cells,IB_ranks,IB_ranks_table,IB_nodes,amr)
    return solid_numbers,IB_ranks,IB_nodes,IB_buffer
end
function initialize_immersed_boundary!(amr::AMR)
    solid_cells = initialize_solidcells(amr)
    solid_numbers,IB_ranks,IB_nodes,IB_buffer = initialize_downwind_transport(solid_cells,amr)
    ib = amr.field.immersed_boundary
    ib.solid_cells = solid_cells
    ib.solid_numbers = solid_numbers
    ib.IB_nodes = IB_nodes
    ib.IB_ranks = IB_ranks
    ib.IB_buffer = IB_buffer
    initialize_solid_neighbor!(amr)
    enc_ghost_target!(amr)
    initialize_upwind2nd!(amr)
    initialize_target_neighbor!(amr)
end
function update_mirror_enc!(ps4est,amr::AMR{DIM,NDF}) where{DIM,NDF}
    GC.@preserve ps4est amr begin
        mirror_data_pointers = amr.ghost.ghost_exchange.mirror_data_pointers
        vs_num = amr.global_data.status.max_vs_num
        gp = PointerWrapper(amr.global_data.forest.ghost)
        pp = PointerWrapper(ps4est)
        for i in eachindex(mirror_data_pointers)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data,InsideSolidData)&&continue
            ap = Base.unsafe_wrap(
                Vector{Cdouble},
                mirror_data_pointers[i],
                3 * DIM + 4 + NDF * vs_num,
            )
            ap[3*DIM+4] = ps_data.bound_enc
        end
    end
end
function bound_enc_exchange!(p4est,amr::AMR{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_enc!(p4est, amr)
    amr_exchange!(
        p4est,
        amr.global_data.forest.ghost,
        amr.ghost.ghost_exchange.ghost_datas,
        amr.ghost.ghost_exchange.mirror_data_pointers,
        size_Ghost_Data(amr.global_data.status.max_vs_num,DIM,NDF),
    )
end
function enc_ghost_target!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    ps4est = amr.global_data.forest.p4est
    update_mirror_enc!(ps4est, amr)
    bound_enc_exchange!(ps4est, amr)
    ghost_wrap = amr.ghost.ghost_wrap
    ghost_datas = amr.ghost.ghost_exchange.ghost_datas
    global_vs_num = amr.global_data.status.max_vs_num
    for i in eachindex(ghost_wrap)
        p = ghost_data_ptr(
            size_Ghost_Data(global_vs_num,DIM,NDF),
            ghost_datas,
            i - 1,
        )
        ds = Base.unsafe_wrap(Vector{Cdouble}, p, DIM)
        ds[1]==EPS&&continue
        offset = 3*DIM+3
        bound_enc = Int(Base.unsafe_load(p + offset * sizeof(Cdouble)))
        ghost_wrap[i].bound_enc = bound_enc
    end
end
function IB_nodes_communicate(solid_cells::Vector{T},IB_ranks::Vector,IB_ranks_table::Vector) where{T<:SolidCells}
    rNumbers,s_vs_nums,r_vs_nums = IB_nodes_pre_communicate(solid_cells,IB_ranks_table) # Communicate Numbers and vs_nums for buffer allocation and topological info quadids.
    IB_buffer = IB_nodes_data_communicate(IB_ranks,s_vs_nums,r_vs_nums)
    MPI.Barrier(MPI.COMM_WORLD)
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
    MPI.Waitall(reqs)
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
    MPI.Waitall(reqs)
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
    MPI.Waitall(reqs)
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
function IB_nodes_data_wrap!(ps4est::P_pxest_t,rNumbers::Vector{Vector{Int}},r_vs_nums::Vector{Vector{Int}},rdata::Vector{IBTransferData},solid_cells::Vector{SolidCells{DIM,NDF}},IB_nodes::Vector,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
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
        for j in eachindex(IB_nodes)
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
                    push!(IB_nodes[j][local_id],ghost_node)
                end
                ps_offset+=1;vs_offset+=vs_num
            end
        end 
    end
    return ghost_nodes
end
function initialize_target_neighbor!(amr)
    ib = amr.field.immersed_boundary
    solid_cells = ib.solid_cells;ib_nodes = ib.IB_nodes
    for i in eachindex(solid_cells)
        for j in eachindex(solid_cells[i].ps_datas)
            initialize_target_neighbor!(solid_cells[i].ps_datas[j],ib_nodes[i][j],amr)
        end
    end
    update_solid_cell!(amr)
end
function initialize_target_neighbor!(ps_data,ib_nodes::Vector{AbstractIBNode},amr::AMR{DIM,NDF}) where{DIM,NDF}
    target_dirs = findall(x->!(isnothing(x[1])||isa(x[1],AbstractInsideSolidData))&&x[1].bound_enc>0,ps_data.neighbor.data[1:2*DIM])
    ib = amr.global_data.config.IB[-ps_data.bound_enc]
    for i in target_dirs
        target_cell = ps_data.neighbor.data[i][1]
        boundary_point,normal = calc_intersect(target_cell.midpoint,ps_data.midpoint,ib)
        downwind_dir,templates,Ainv = initialize_templates(target_cell,ib_nodes,boundary_point,normal,i)
        ps_data.neighbor.data[i][1] = TargetNeighbor{DIM,NDF}(i,boundary_point,normal,downwind_dir,templates,Ainv,target_cell)
    end
end
function initialize_templates(target_cell::AbstractPsData{2,NDF},ib_nodes::Vector{AbstractIBNode},boundary_point,normal::Vector{Float64},faceid::Int)where{NDF}
    DIM = 2
    dir = get_dir(faceid);rot = get_rot(faceid)
    face_normal = zeros(DIM);face_normal[dir] = -rot
    Rot = @SMatrix [0. -1.;1. 0.]
    rdir = sign(normal[1]*face_normal[2]-normal[2]*face_normal[1])
    downwind_dir = Array(rdir*Rot*(face_normal+normal))
    downwind_dir./=norm(downwind_dir)
    weights = zeros(length(ib_nodes))
    for i in eachindex(ib_nodes)
        r = boundary_point-ib_nodes[i].midpoint
        w = dot(downwind_dir,r)
        if ib_nodes[i].midpoint==target_cell.midpoint
            weights[i] = -Inf
        elseif  w>0
            weights[i] = w/norm(r)^2
        else
            weights[i] = w
        end
    end
    indices = sortperm(weights;rev=true)
    upwind_nodes = ib_nodes[indices]
    # distance = [norm(x.midpoint-boundary_point) for x in ib_nodes]
    templates,A = make_A_2D!(upwind_nodes)
    # if target_cell.midpoint[1]>0.105-EPS&&target_cell.midpoint[1]<0.105+EPS&&target_cell.midpoint[2]>0.995-EPS&&target_cell.midpoint[2]<0.995+EPS
    #     @show A
    # end
    # if target_cell.midpoint[1]>0.105-EPS&&target_cell.midpoint[1]<0.105+EPS&&target_cell.midpoint[2]>-0.995-EPS&&target_cell.midpoint[2]<-0.995+EPS
    #     @show A
    # end
    return downwind_dir,templates,inv(A)
end
function make_basis(x)
    if length(x)==2
        return [x[1]^2,x[2]^2,x[1]*x[2],x[1],x[2],1.0]
    else
        return [x[1]^2,x[2]^2,x[3]^2,x[1]*x[2],x[1]*x[3],x[2]*x[3],x[1],x[2],x[3],1.0]
    end
end
function make_slope_basis(x,dir)
    if length(x)==2
        if dir==1
            return [2.0*x[1],0.,x[2],1.,0.,0.]
        else
            return [0.,2.0*x[2],x[1],0.,1.,0.]
        end
    else
        if dir==1
            return [2.0*x[1],0.,0.,x[2],x[3],0.,1.,0.,0.,0.]
        elseif dir==2
            return [0.,2.0*x[2],0.,x[1],0.,x[3],0.,1.,0.,0.]
        else
            return [0.,0.,2.0*x[3],0.,x[1],x[2],0.,0.,1.,0.]
        end
    end
end
function make_A_2D!(nodes::Vector{AbstractIBNode})
    k = 0
    templates = Vector{AbstractIBNode}(undef,6)
    selected = nothing
    for i in eachindex(nodes)
        if k == 0
            selected = make_basis(nodes[i].midpoint)
            k = 1
            templates[k] = nodes[i]
        else
            temp = hcat(selected,make_basis(nodes[i].midpoint))
            r = rank(temp)
            if r > k
                selected = temp
                k = r
                templates[k] = nodes[i]
            end
        end
        if k == 6
            break
        end
    end
    if k<6
        @show k length(nodes)
        throw(`Larger search radius is required!`)
    end
    return templates,permutedims(selected)
end
function update_solid_cell!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    IB_data_exchange!(amr)
    scs = amr.field.immersed_boundary.solid_cells
    for solid_cells in scs 
        for solid_cell in solid_cells.ps_datas
            target_dirs = findall(x->isa(x[1],TargetNeighbor),solid_cell.neighbor.data)
            vs_data = solid_cell.vs_data
            b = make_basis(solid_cell.midpoint)
            C = Array{Float64}(undef,length(b),vs_data.vs_num,NDF)
            if length(target_dirs)>1
                df = vs_data.df;sdf = vs_data.sdf
                df.=0.
                weights = Matrix{Float64}(undef,vs_data.vs_num,length(target_dirs))
                for j in axes(weights,2)
                    downwind_dir = solid_cell.neighbor.data[target_dirs[j]][1].downwind_dir
                    for i in axes(weights,1)
                        @views weights[i,j] = max(0.,dot(downwind_dir,vs_data.midpoint[i,:]))
                    end
                end
                weight_sum = sum(weights;dims=2)
                for i in axes(weights,1)
                    @views weights[i,:] .= weight_sum[i]==0. ? 0. : weights[i,:]./weight_sum[i]
                end
                for ti in eachindex(target_dirs)
                    target_neighbor = solid_cell.neighbor.data[target_dirs[ti]][1]
                    dir = get_dir(target_neighbor.faceid)
                    for i in eachindex(target_neighbor.templates)
                        vs_data_t = target_neighbor.templates[i].vs_data                    
                        @views vs_project!(vs_data_t.df,vs_data_t.level,C[i,:,:],vs_data.level,vs_data)
                    end
                    Ainv = target_neighbor.Ainv
                    sb = make_slope_basis(target_neighbor.target_cell.midpoint,get_dir(target_neighbor.faceid))
                    for i in axes(df,1)
                        weights[i,ti]==0&&continue
                        for j in axes(df,2)
                            c = @views Ainv*C[:,i,j]
                            df[i,j]+=dot(b,c)*weights[i,ti]
                            sdf[i,j,dir]=dot(sb,c)
                        end
                    end
                end
            else
                target_neighbor = solid_cell.neighbor.data[target_dirs[1]][1]
                dir = get_dir(target_neighbor.faceid)
                for i in eachindex(target_neighbor.templates)
                    vs_data_t = target_neighbor.templates[i].vs_data                    
                    @views vs_project!(vs_data_t.df,vs_data_t.level,C[i,:,:],vs_data.level,vs_data)
                end
                df = vs_data.df;sdf = vs_data.sdf
                Ainv = target_neighbor.Ainv
                sb = make_slope_basis(target_neighbor.target_cell.midpoint,get_dir(target_neighbor.faceid))
                for i in axes(df,1)
                    for j in axes(df,2)
                        c = @views Ainv*C[:,i,j]
                        df[i,j]=dot(b,c)
                        sdf[i,j,dir]=dot(sb,c)
                    end
                end
            end
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
function initialize_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_dirs = findall(x->!isnothing(x[1])&&x[1].bound_enc<0,ps_data.neighbor.data[1:2*DIM])
    vs_data = ps_data.vs_data
    for i in solid_dirs
        solid_cell = ps_data.neighbor.data[i][1]
        if ps_data.bound_enc!=0
            ps_data.bound_enc!=-solid_cell.bound_enc&&throw(`Different immersed boundaries correspond to the same fluid cell!`)
        end
        ps_data.bound_enc=-solid_cell.bound_enc
        ib = amr.global_data.config.IB[-solid_cell.bound_enc]
        aux_point,normal = calc_intersect(ps_data.midpoint,solid_cell.midpoint,ib)
        svsdata = VS_Data{DIM,NDF}(
            vs_data.vs_num,
            vs_data.level,
            vs_data.weight,
            vs_data.midpoint,
            zeros(vs_data.vs_num,NDF),
            zeros(vs_data.vs_num,NDF,DIM),
            Matrix{Float64}(undef,0,0)
        )
        cvc = initialize_cutted_velocity_cell(normal,svsdata,amr) # heavy overhead
        ps_data.neighbor.data[i][1] = SolidNeighbor{DIM,NDF}(
            solid_cell.bound_enc,i,aux_point,normal,solid_cell,solid_cell.midpoint,solid_cell.ds,zeros(DIM+2),zeros(DIM+2,DIM),zeros(vs_data.vs_num,NDF),cvc,svsdata
        )
    end
    return nothing
end
function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        update_solid_neighbor!(amr.global_data.config.solver.flux,ps_data,ps_data.neighbor.data[i][1],amr)
    end
end
function update_solid_neighbor!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    update_upwind2nd!(amr)
    for tree in amr.field.trees.data
        for ps_data in tree
            # if !isa(ps_data,InsideSolidData)&&ps_data.midpoint[1]>0.105-EPS&&ps_data.midpoint[1]<0.105+EPS&&ps_data.midpoint[2]>0.995-EPS&&ps_data.midpoint[2]<0.995+EPS
            #     @show ps_data.midpoint ps_data.bound_enc
            # end
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
            update_solid_neighbor!(ps_data,amr)
        end
    end
end
function update_solid_neighbor!(::AbstractFluxType,ps_data::PS_Data{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF},amr::AMR) where{DIM,NDF}
    global_data = amr.global_data;ib = global_data.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid);rot = get_rot(faceid)
    dx = solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    opp = faceid%2==0 ? faceid-1 : faceid+1
    upwind1 = ps_data.neighbor.data[opp][1];vs_data1 = upwind1.vs_data
    df1 = vs_project(vs_data1.df,vs_data1.level,vs_data.level,vs_data)
    Θ = heaviside.(vn) # Solid to gas (emission) is 1.
    Θ_grid = @views heaviside.(rot*vs_data.midpoint[:,dir]) # Solid to gas (downwind) is 1.
    construct_incident_df!(aux_df,ps_data,solid_neighbor;Θ,Θ_grid)
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
    Mu_L,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
    ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
    aux_prim = IB_prim(ib,aux_point,ρw)
    F = aux_prim[1]*M
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            aux_df[i,:] .= F[i,:]
        end
    end
    cvc_correction!(aux_df,F,solid_neighbor,amr)
    x_ip,df_ip = reconstruct_image_point(ps_data,solid_neighbor;df1) # Reconstruct 2nd order df with upwind2nd_df at the image point. x_ip is dimensionless.
    x_b = aux_point[dir];x = ps_data.midpoint[dir]
    A = @SMatrix [(x-dx)^2 x-dx 1.;x_ip^2 x_ip 1.;x_b^2 x_b 1.]
    Ainv = inv(A)
    reconstruct_solid_neighbor!(Ainv,df1,df_ip,aux_df,
        solid_neighbor.vs_data.df;Θ,Θ_grid,ps_data,solid_neighbor,dir,dx)
    reconstruct_solid_slope!(ps_data,solid_neighbor;dir,dx,Θ,Θ_grid)
    update_target_sdf!(ps_data,solid_neighbor;Θ_grid)
    # if ps_data.midpoint[1]>0.105-EPS&&ps_data.midpoint[1]<0.105+EPS&&ps_data.midpoint[2]>0.995-EPS&&ps_data.midpoint[2]<0.995+EPS
    #     @show `up` maximum(abs.(vs_data.sdf[:,:,dir]))
    # end
    # if ps_data.midpoint[1]>0.105-EPS&&ps_data.midpoint[1]<0.105+EPS&&ps_data.midpoint[2]>-0.995-EPS&&ps_data.midpoint[2]<-0.995+EPS
    #     @show `down` maximum(abs.(vs_data.sdf[:,:,dir]))
    # end
    solid_neighbor.w = calc_w0(vs_data.midpoint,solid_neighbor.vs_data.df,vs_data.weight,global_data)
    solid_neighbor.sw[:,dir] .= (solid_neighbor.w-ps_data.w)./(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])
end
function reconstruct_image_point(ps_data::PS_Data{2,NDF},solid_neighbor;kwargs...) where{NDF}
    df1 = kwargs[:df1]
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid)
    vs_data = ps_data.vs_data
    x = ps_data.midpoint[dir]
    dx = solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]
    # x_ip = (solid_neighbor.aux_point[dir]-0.5*dx-ps_data.midpoint[dir])/dx # The origin is set at target cell's center.
    df2 = solid_neighbor.upwind2nd_df;df = vs_data.df
    # df_ip = @. (0.5*df-df1+0.5*df2)*x_ip^2+(1.5*df-2.0*df1+0.5*df2)*x_ip+df
    x_ip = solid_neighbor.aux_point[dir]-0.5*dx
    A = @SMatrix [(x-2*dx)^2 (x-2*dx) 1.0;(x-dx)^2 (x-dx) 1.0;x^2 x 1.0]
    Ainv = inv(A)
    v = @SVector [x_ip^2,x_ip,1.0]
    df_ip = Matrix{Float64}(undef,vs_data.vs_num,NDF)
    b = zeros(3)
    for j in axes(df,2)
        for i in axes(df,1)
            b[1] = df2[i,j];b[2] = df1[i,j];b[3] = df[i,j]
            c = Ainv*b
            df_ip[i,j] = dot(v,c)
        end
    end
    return x_ip,df_ip
end
function reconstruct_solid_neighbor!(Ainv,df1,df_ip,aux_df,solid_df;ps_data::PS_Data{DIM,NDF},kwargs...) where{DIM,NDF}
    dir = kwargs[:dir];dx = kwargs[:dx];Θ = kwargs[:Θ];Θ_grid = kwargs[:Θ_grid]
    x = ps_data.midpoint[dir]
    basis = @SVector [(x+dx)^2,x+dx,1.0]
    solid_neighbor = kwargs[:solid_neighbor]
    vs_data = ps_data.vs_data;svs_data = solid_neighbor.solid_cell.vs_data
    sdf_2nd = @views solid_neighbor.vs_data.sdf[:,:,dir]
    level = vs_data.level;nlevel = svs_data.level
    b = Vector{Float64}(undef,3)
    flag = 0.;index = 1
    for i in axes(solid_df,1)
        if Θ_grid[i]==1.
            if Θ[i]==1. # solid part
                for j in axes(solid_df,2)
                    b[1] = df1[i,j];b[2] = df_ip[i,j];b[3] = aux_df[i,j]
                    c = Ainv*b
                    solid_df[i,j] = dot(c,basis) # Dimensionless position of the solid_neighbor is 1.0, which corresponds to the basis [1.,1.,1.]
                    sdf_2nd[i,j] = 2*c[1]*x+c[2] # Second order derivatives at the target cell for the update of solid_neighbor's sdf
                end
                if level[i]==nlevel[index]
                    index+=1
                elseif level[i]<nlevel[index]
                    while flag != 1.0
                        flag += 1/2^(DIM*(nlevel[index]-level[i]))
                        index += 1
                    end
                    flag = 0.
                else
                    flag += 1/2^(DIM*(level[i]-nlevel[index]))
                    if flag == 1.
                        index += 1
                        flag = 0.
                    end
                end
            else # downwind part
                if level[i]==nlevel[index]
                    for j in axes(solid_df,2)
                        solid_df[i,j] = svs_data.df[index,j]
                    end
                    index+=1
                elseif level[i]<nlevel[index]
                    solid_df[i,:] .= 0.
                    while flag != 1.0
                        for j in axes(solid_df,2)
                            solid_df[i,j] += svs_data.df[index,j]/2^(DIM*(nlevel[index]-level[i]))
                        end
                        flag += 1/2^(DIM*(nlevel[index]-level[i]))
                        index += 1
                    end
                    flag = 0.
                else
                    flag += 1/2^(DIM*(level[i]-nlevel[index]))
                    for j in axes(solid_df,2)
                        solid_df[i,j] = svs_data.df[index,j]
                    end
                    if flag == 1.
                        index += 1
                        flag = 0.
                    end
                end
            end
        else
            if level[i]==nlevel[index]
                index+=1
            elseif level[i]<nlevel[index]
                while flag != 1.0
                    flag += 1/2^(DIM*(nlevel[index]-level[i]))
                    index += 1
                end
                flag = 0.
            else
                flag += 1/2^(DIM*(level[i]-nlevel[index]))
                if flag == 1.
                    index += 1
                    flag = 0.
                end
            end
        end
    end
end
function reconstruct_solid_slope!(ps_data::PS_Data{DIM,NDF},solid_neighbor;kwargs...) where{DIM,NDF}
    dir = kwargs[:dir];dx = kwargs[:dx];Θ = kwargs[:Θ];Θ_grid = kwargs[:Θ_grid]
    scsdf = @views solid_neighbor.solid_cell.vs_data.sdf[:,:,dir] # sdf_2nd for downwind part
    snsdf = @views solid_neighbor.vs_data.sdf[:,:,dir]# sdf_2nd for solid part
    df_ex = solid_neighbor.vs_data.df
    df = ps_data.vs_data.df
    level = ps_data.vs_data.level;nlevel = solid_neighbor.solid_cell.vs_data.level
    flag = 0.;index = 1
    for i in axes(snsdf,1)
        if Θ_grid[i]==1.
            if Θ[i]==1. # solid part
                for j in axes(snsdf,2)
                    snsdf[i,j] = -2.0*snsdf[i,j]-3.0*(df[i,j]-df_ex[i,j])/dx
                end
                if level[i]==nlevel[index]
                    index+=1
                elseif level[i]<nlevel[index]
                    while flag != 1.0
                        flag += 1/2^(DIM*(nlevel[index]-level[i]))
                        index += 1
                    end
                    flag = 0.
                else
                    flag += 1/2^(DIM*(level[i]-nlevel[index]))
                    if flag == 1.
                        index += 1
                        flag = 0.
                    end
                end
            else # downwind part
                if level[i]==nlevel[index]
                    for j in axes(snsdf,2)
                        snsdf[i,j] = -2.0*scsdf[index,j]-3.0*(df[i,j]-df_ex[i,j])/dx
                    end
                    index+=1
                elseif level[i]<nlevel[index]
                    for j in axes(snsdf,2)
                        snsdf[i,j] = -3.0*(df[i,j]-df_ex[i,j])/dx
                    end
                    while flag != 1.0
                        for j in axes(snsdf,2)
                            snsdf[i,j]-=2.0*scsdf[index,j]/2^(DIM*(nlevel[index]-level[i]))
                        end
                        flag += 1/2^(DIM*(nlevel[index]-level[i]))
                        index += 1
                    end
                    flag = 0.
                else
                    for j in axes(snsdf,2)
                        snsdf[i,j] = -2.0*scsdf[index,j]-3.0*(df[i,j]-df_ex[i,j])/dx
                    end
                    flag += 1/2^(DIM*(level[i]-nlevel[index]))
                    if flag == 1.
                        index += 1
                        flag = 0.
                    end
                end
            end
        else
            if level[i]==nlevel[index]
                index+=1
            elseif level[i]<nlevel[index]
                while flag != 1.0
                    flag += 1/2^(DIM*(nlevel[index]-level[i]))
                    index += 1
                end
                flag = 0.
            else
                flag += 1/2^(DIM*(level[i]-nlevel[index]))
                if flag == 1.
                    index += 1
                    flag = 0.
                end
            end
        end
    end
end
function update_target_sdf!(ps_data,solid_neighbor;kwargs...)
    Θ_grid = kwargs[:Θ_grid]
    dir = get_dir(solid_neighbor.faceid)
    df = ps_data.vs_data.df;sndf = solid_neighbor.vs_data.df
    sdf = @views ps_data.vs_data.sdf[:,:,dir]
    dx = solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]
    for j in axes(sdf,2)
        for i in axes(sdf,1)
            if Θ_grid[i]==1.0
                sdf[i,j] = (sndf[i,j]-df[i,j])/dx
            end
        end
    end
end
function construct_incident_df!(aux_df,ps_data::PS_Data{DIM,NDF},solid_neighbor;kwargs...) where{DIM,NDF}
    Θ = kwargs[:Θ];Θ_grid = kwargs[:Θ_grid]
    dir = get_dir(solid_neighbor.faceid)
    dx = solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]
    dxL = solid_neighbor.aux_point[dir]-ps_data.midpoint[dir]
    vs_data = ps_data.vs_data
    scvs_data = solid_neighbor.solid_cell.vs_data
    df = vs_data.df
    level = vs_data.level;nlevel = scvs_data.level
    flag = 0.;index = 1
    sdf = @views vs_data.sdf[:,:,dir];scsdf = @views scvs_data.sdf[:,:,dir]
    for i in axes(df,1)
        if Θ_grid[i]==1.0&&Θ[i]==0. # downwind part
            if level[i]==nlevel[index]
                for j in axes(df,2)
                    aux_df[i,j] = df[i,j]+(scvs_data.df[index,j]-df[i,j])/dx*dxL
                end
                index+=1
            elseif level[i]<nlevel[index]
                while flag != 1.0
                    for j in axes(df,2)
                        aux_df[i,j] += (df[i,j]+scsdf[index,j]*dxL)/2^(DIM*(nlevel[index]-level[i]))
                    end
                    flag += 1/2^(DIM*(nlevel[index]-level[i]))
                    index += 1
                end
                flag = 0.
            else
                for j in axes(df,2)
                    aux_df[i,j] = df[i,j]+scsdf[index,j]*dxL
                end
                flag += 1/2^(DIM*(level[i]-nlevel[index]))
                if flag == 1.
                    index += 1
                    flag = 0.
                end
            end
        else
            for j in axes(df,2)
                aux_df[i,j] = df[i,j]+sdf[i,j]*dxL
            end
            if level[i]==nlevel[index]
                index+=1
            elseif level[i]<nlevel[index]
                while flag != 1.0
                    flag += 1/2^(DIM*(nlevel[index]-level[i]))
                    index += 1
                end
                flag = 0.
            else
                flag += 1/2^(DIM*(level[i]-nlevel[index]))
                if flag == 1.
                    index += 1
                    flag = 0.
                end
            end
        end
    end
end