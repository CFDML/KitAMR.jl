# function boundary_interpolation_points(amr::AMR,nodes::AbstractMatrix,r::Real) # Search fluid points inside the detect radius r
#     interp_points = Vector{Float64}
#     datas = amr.field.trees.data
#     for i in eachindex(data)
#         for j in eachindex(data[i])
#             ps_data = data[i][j]

#         end
#     end
# end

function solid_flag(boundary::Circle,midpoint::AbstractVector) # Does midpoint locate at solid?
    return xor(norm(midpoint.-boundary.center)>boundary.radius,boundary.solid)
end
function solid_flag(::Domain{Maxwellian},::AbstractVector) # Does midpoint locate at solid?
    return false
end
function solid_cell_flag(boundary::Circle,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data,inside::Bool) # Screen ghost cells, that is those inside solid domain and immediately adjacent the boundary.
    (boundary_flag(boundary,midpoint,ds,global_data) && xor(boundary.solid,!inside)) && return true
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
function calc_intersect_point(boundary::Circle,midpoint::AbstractVector{Float64})
    boundary.radius/norm(midpoint-boundary.center)*(midpoint-boundary.center)+boundary.center
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
function broadcast_boundary_midpoints!(boundary_points::Vector{Vector{Vector{Float64}}},::Global_Data{DIM})where{DIM}
    Numbers = pre_broadcast_boundary_points(boundary_points)
    rbuffer = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Float64}}(undef,length(boundary_points)) # boundaries{points}
    for i in eachindex(boundary_points)
        buffer = Vector{Float64}(undef,DIM*length(boundary_points[i]))
        for j in eachindex(boundary_points[i])
            buffer[DIM*(i-1)+1:DIM*i] .= boundary_points[i][j]
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
    for i in eachindex(boundary_points)
        boundary_points_global[i] = Vector{Float64}[]
        for j in eachindex(Numbers)
            if j-1 == MPI.Comm_rank(MPI.COMM_WORLD)
                append!(boundary_points_global[i],boundary_points[i])
            else
                for k in 1:Int(length(rbuffer[i][j])/DIM)
                    push!(boundary_points_global[i],rbuffer[i][j][(DIM*(k-1)+1):DIM*k]) 
                end
            end
        end
    end
    # for i in eachindex(solid_cells)
    #     buffer = Vector{Float64}(undef,DIM*length(solid_cells[i].ps_datas))
    #     for j in eachindex(solid_cells[i].ps_datas)
    #         buffer[DIM*(i-1)+1:DIM*i] .= aux_points[i].midpoints[j]
    #     end
    #     sbuffer[i] = buffer
    # end
    # for i in eachindex(solid_cells)
    #     for j in eachindex(Numbers)
    #         j-1==MPI.Comm_rank(MPI.COMM_WORLD) && (rbuffer[i][j] = sbuffer[i])
    #     end
    # end
    # for i in eachindex(solid_cells)
    #     for j in eachindex(Numbers)
    #         MPI.Bcast!(rbuffer[i][j],j-1,MPI.COMM_WORLD)
    #     end
    # end
    # MPI.Barrier(MPI.COMM_WORLD)
    # for i in eachindex(solid_cells)
    #     for j in eachindex(Numbers)
    #         j-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
    #         for k in 1:length(rbuffer[i][j])/DIM
    #            push!(aux_points.midpoints,rbuffer[i][j][DIM*(k-1)+1:DIM*k])
    #            push!(aux_points.mpi_rank,j-1)
    #         end
    #     end
    # end
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
end

function IB_flag(boundary::Circle,aux_point::AbstractVector,midpoint::AbstractVector,ds::AbstractVector)
    DIM = length(aux_point)
    r = boundary.search_radius
    for i in eachindex(2^DIM)
        norm(midpoint.+0.5*ds.*RMT[DIM][i].-aux_point)<r && return true # any corner cross boundary?
    end
    return false
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

function search_IB!(IB_nodes::Vector{Vector{Vector{PS_Data{DIM,NDF}}}},aux_points::Vector,trees::PS_Trees{DIM,NDF},global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    ps_datas = trees.data
    boundaries = global_data.config.IB
    for i in eachindex(ps_datas)
        for j in eachindex(ps_datas[i])
            ps_data = ps_datas[i][j]
            isa(ps_data,InsideSolidData) && continue
            for k in eachindex(boundaries)
                for l in eachindex(aux_points[k])
                    IB_flag(boundaries[k],aux_points[k][l],ps_data.midpoint,ps_data.ds)&&push!(IB_nodes[k][l],ps_data)
                end
            end
        end
    end
end

function IB_data_exchange!(boundary::Boundary)
    IB_ranks_table = boundary.IB_ranks_table
    sdata = boundary.IB_buffer.sdata
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    s_data_buffer = Vector{Matrix{Cdouble}}(undef,length(sdata))
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(IB_ranks_table)
        i-1==rank&&continue
        isempty(IB_ranks_table[i])&&continue
        for j in eachindex(IB_ranks_table[i])
            vnums[i]+=IB_ranks_table[i][j].vs_data.vs_num
        end
        ps_num = length(IB_ranks_table[i])
        s_data_buffer[i] = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums[i],NDF))
        offset = 0
        for j in eachindex(IB_ranks_table[i])
            vs_num = IB_ranks_table[i][j].vs_data.vs_num
            s_data_buffer[i][offset+1:offset+vs_num,:] .= IB_ranks_table[i][j].vs_data.df
            offset+=vs_num
        end
        sreq = MPI.Isend(
            s_data_buffer[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    rdata = boundary.IB_buffer.rdata
    r_data_buffer = Vector{Matrix{Cdouble}}(undef,length(sdata))
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    for i in eachindex(r_data_buffer)
        i-1==rank&&continue
        !isdefined(r_vs_nums,i)&&continue
        vnums = sum(r_vs_nums[i])
        r_data_buffer[i] = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(rdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums,NDF))
        rreq = MPI.Irecv!(
                r_data_buffer[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
end
function IB_vs_nums_exchange!(boundary::Boundary)
    IB_ranks_table = boundary.IB_ranks_table
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    sbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{solid_cells{vs_data{}}}
    for i in eachindex(sbuffer)
        i-1 == rank&&continue
        sbuffer[i] = Vector{Int}(undef,length(IB_ranks_table[i]))
        for j in eachindex(IB_ranks_table[i])
            sbuffer[i][j] = IB_ranks_table[i][j].vs_data.vs_num
        end
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
        !isdefined(r_vs_nums,i) && continue
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
function IB_structure_exchange!(boundary::Boundary)
    IB_ranks_table = boundary.IB_ranks_table
    sdata = boundary.IB_buffer.sdata
    for i in eachindex(sdata)
        sc_free(-1,sdata[i])
    end
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(IB_ranks_table)
        i-1==rank&&continue
        isempty(IB_ranks_table[i])&&continue
        for j in eachindex(IB_ranks_table[i])
            vnums[i]+=IB_ranks_table[i][j].vs_data.vs_num
        end
        ps_num = length(IB_ranks_table[i])
        sdata[i] = sc_malloc(-1,vnums[i]*(sizeof(Cdouble)*NDF+sizeof(Int8))+ps_num*DIM*sizeof(Cdouble))
        midpoint_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]),(ps_num,DIM))
        df_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums[i],NDF))
        level_buffer = unsafe_wrap(Vector{Int8},Ptr{Int8}(sdata[i]+ps_num*DIM*sizeof(Cdouble)+vnums[i]*sizeof(Cdouble)*NDF),vnums[i])
        offset = 0
        for j in eachindex(IB_ranks_table[i])
            midpoint_buffer[j,:] .= IB_ranks_table[i][j].midpoint
        end
        for j in eachindex(IB_ranks_table[i])
            vs_num = IB_ranks_table[i][j].vs_data.vs_num
            df_buffer[offset+1:offset+vs_num,:] .= IB_ranks_table[i][j].vs_data.df
            level_buffer[offset+1:offset+vs_num] .= IB_ranks_table[i][j].vs_data.level
            offset+=vs_num
        end
        sreq = MPI.Isend(
            sdata[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + rank,
        )
        push!(reqs, sreq)
    end
    rdata = boundary.IB_buffer.rdata
    for i in eachindex(rdata)
        sc_free(-1,rdata[i])
    end
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    for i in eachindex(rdata)
        i-1 == rank&&continue
        !isdefined(r_vs_nums,i) && continue
        ps_nums = length(r_vs_nums[i])
        rdata[i] = sc_malloc(-1,sum(r_vs_nums[i])*(sizeof(Cdouble)*NDF+sizeof(Int8))+DIM*ps_nums*sizeof(Cdouble))
        rreq = MPI.Irecv!(
                rdata[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
end
function IB_wrap_update!(boundary::Boundary)
    rdata = boundary.IB_buffer.rdata
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(rdata)
        i-1==rank&&continue
        !isdefined(r_vs_nums,i)&&continue
        ps_offset = 1;vs_offset = 0
        vnums = sum(r_vs_nums[i])
        midpoints = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(rdata[i]),(length(r_vs_nums[i]),DIM))
        df_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(rdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums,NDF))
        level_buffer = unsafe_wrap(Vector{Int8},Ptr{Int8}(rdata[i]+ps_num*DIM*sizeof(Cdouble)+vnums*sizeof(Cdouble)*NDF),vnums)
        for j in eachindex(IB_cells)
            for k in eachindex(IB_cells[j].IB_nodes)
                for l in eachindex(IB_cells[j].IB_nodes[k])
                    !isa(IB_cells[j].IB_nodes[k][l],GhostIBNode) && continue
                    vs_num = r_vs_nums[i][ps_offset]
                    node = IB_cells[j].IB_nodes[k][l]
                    node.midpoint = @view(midpoints[ps_offset,:])
                    node.level = @view(level_buffer[vs_offset+1:vs_offset+vs_num])
                    node.df = @view(df_buffer[vs_offset+1:vs_offset+vs_num,:])
                    ps_offset+=1;vs_offset+=vs_num
                end
            end
        end
    end
end
function IB_structure_update!(boundary::Boundary)
    IB_vs_nums_exchange!(boundary)
    IB_structure_exchange!(boundary)
    IB_wrap_update!(boundary)
end
function IB_quadid_update!(boundary::Boundary) # Before partition, global quadids of solid_cells and IB_nodes need to be updated
    broadcast_quadid!(boundary.solid_cells)
end
function IB_update!(boundary::Boundary) # After partition, IB need to be reinitialized
    
end