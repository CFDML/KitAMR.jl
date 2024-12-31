get_bc(bc::AbstractVector) = bc
function solid_cell_index_encoder!(solid_cell_index::Vector{Int},now_index::Int)
    # if solid_cell_index>2^120
    #     @error `Too many solid cells are sharing the same IB node!`
    # end
    # if now_index>2^20
    #     @error `Too many solid cells that a larger index type is needed!`
    # end
    # if solid_cell_index!=0
    #     solid_cell_index = solid_cell_index<<20 + now_index
    # else
    #     solid_cell_index = now_index
    # end
    id = findfirst(x->x==0,solid_cell_index)
    isnothing(id) && (@error `A larger SOLID_CELL_ID_NUM is needed!`)
    solid_cell_index[id]=now_index
end
function solid_cell_index_decoder(solid_cell_index::Vector{Int})
    # ids = Int[]
    # for _ in 1:6
    #     a = solid_cell_index%2^20
    #     a==0&&break
    #     push!(ids,a)
    #     solid_cell_index = solid_cell_index>>20
    # end
    # return ids
    ids = findall(x->x!=0,solid_cell_index)
    @show solid_cell_index[ids]
    return solid_cell_index[ids]
end
function solid_flag(boundary::Circle,midpoint::AbstractVector) # Does midpoint locate at solid?
    return xor(norm(midpoint.-boundary.center)>boundary.radius,boundary.solid)
end
function solid_flag(::Domain{Maxwellian},::AbstractVector) # Does midpoint locate at solid?
    return false
end
function solid_cell_flag(boundary::Circle,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data,inside::Bool) # Ghost nodes, those are inside solid domain and immediately adjacent the boundary.
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

function IB_flag(boundary::Circle,aux_point::AbstractVector,midpoint::AbstractVector,::AbstractVector)
    # DIM = length(aux_point)
    r = boundary.search_radius
    #=
    for i in eachindex(2^DIM)
	    corner = midpoint.+0.5*ds.*RMT[DIM][i]
	    (norm(corner.-aux_point)<r||norm(corner-boundary.center)<boundary.radius) && return true # any corner cross boundary?
    end
    return false
    =#
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
function search_IB!(IB_nodes::Vector{Vector{Vector{PS_Data{DIM,NDF}}}},aux_points::Vector,trees::PS_Trees{DIM,NDF},global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
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
                        solid_cell_index_encoder!(ps_data.solid_cell_index,l)
                        push!(IB_nodes[k][l],ps_data)
                    end
                end
            end
        end
    end
end

function update_IB_sdata!(boundary::Boundary)
    sdata = boundary.IB_buffer.sdata
    IB_ranks_table = boundary.IB_ranks_table
    for i in eachindex(IB_ranks_table)
        i-1==MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(IB_ranks_table[i])&&continue
        offset = 0
        for j in eachindex(IB_ranks_table[i])
            sdata[i].prims[j,:] .= IB_ranks_table[i][j].prim
        end
        for j in eachindex(IB_ranks_table[i])
            vs_num = IB_ranks_table[i][j].vs_data.vs_num
            sdata[i].df[offset+1:offset+vs_num,:] .= IB_ranks_table[i][j].vs_data.df
            offset+=vs_num
        end
    end
end
function IB_data_exchange!(amr::AMR)
    boundary = amr.field.boundary
	update_IB_sdata!(boundary)
    IB_ranks_table = boundary.IB_ranks_table
    sdata = boundary.IB_buffer.sdata
    rdata = boundary.IB_buffer.rdata
    r_vs_nums = boundary.IB_buffer.r_vs_nums
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    # s_data_buffer = Vector{Matrix{Cdouble}}(undef,length(sdata))
    # reqs = Vector{MPI.Request}(undef, 0)
    # vnums = zeros(Int,length(IB_ranks_table))
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sdata)
        i-1 == rank&&continue
        isempty(IB_ranks_table[i])&&continue
        # sreq = MPI.Isend(
        #     sdata[i].midpoints,
        #     MPI.COMM_WORLD;
        #     dest = i - 1,
        #     tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        # )
        # push!(reqs, sreq)
        sreq = MPI.Isend(
            sdata[i].prims,
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + rank,
        )
        push!(reqs, sreq)
        # sreq = MPI.Isend(
        #     sdata[i].level,
        #     MPI.COMM_WORLD;
        #     dest = i - 1,
        #     tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        # )
        # push!(reqs, sreq)
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
        # ps_num = length(r_vs_nums[i])
        # vs_nums = sum(r_vs_nums[i])
        # midpoints = Matrix{Cdouble}(undef,ps_num,DIM)
        # prims = Matrix{Cdouble}(undef,ps_num,DIM+2)
        # df = Matrix{Cdouble}(undef,vs_nums,NDF)
        # level = Vector{Int8}(undef,vs_nums)
        # rdata[i] = IBTransferData(midpoints,prims,level,df)
        # rdata[i] = sc_malloc(-1,sum(r_vs_nums[i])*(sizeof(Cdouble)*NDF+sizeof(Int8))+(2*DIM+2)*ps_nums*sizeof(Cdouble))
        # rreq = MPI.Irecv!(
        #         midpoints,
        #         MPI.COMM_WORLD;
        #         source = i - 1,
        #         tag = COMM_DATA_TAG + i - 1,
        #     )
        # push!(reqs, rreq)
        rreq = MPI.Irecv!(
                rdata[i].prims,
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        # rreq = MPI.Irecv!(
        #         level,
        #         MPI.COMM_WORLD;
        #         source = i - 1,
        #         tag = COMM_DATA_TAG + i - 1,
        #     )
        # push!(reqs, rreq)
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
# function IB_slope_exchange!(amr::AMR)
#     boundary = amr.field.boundary
#     IB_ranks_table = boundary.IB_ranks_table
#     sdata = boundary.IB_buffer.sdata
#     rdata = boundary.IB_buffer.rdata
#     r_vs_nums = boundary.IB_buffer.r_vs_nums
#     rank = MPI.Comm_rank(MPI.COMM_WORLD)
#     reqs = Vector{MPI.Request}(undef, 0)
#     for i in eachindex(sdata)
#         i-1 == rank&&continue
#         isempty(IB_ranks_table[i])&&continue
#         sreq = MPI.Isend(
#             sdata[i].sdf,
#             MPI.COMM_WORLD;
#             dest = i - 1,
#             tag = COMM_DATA_TAG + rank,
#         )
#         push!(reqs, sreq)
#     end
#     for i in eachindex(rdata)
#         i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
#         !isassigned(r_vs_nums,i) && continue
#         rreq = MPI.Irecv!(
#                 rdata[i].sdf,
#                 MPI.COMM_WORLD;
#                 source = i - 1,
#                 tag = COMM_DATA_TAG + i - 1,
#             )
#         push!(reqs, rreq)
#     end
#     MPI.Waitall!(reqs)
# end
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
function IB_structure_exchange!(boundary::Boundary)
    IB_ranks_table = boundary.IB_ranks_table
    sdata = boundary.IB_buffer.sdata
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(IB_ranks_table)
        i-1==rank&&continue
        isempty(IB_ranks_table[i])&&continue
        for j in eachindex(IB_ranks_table[i])
            vnums[i]+=IB_ranks_table[i][j].vs_data.vs_num
        end
        ps_num = length(IB_ranks_table[i])
        vs_nums = vnums[i]
        midpoints = Matrix{Cdouble}(undef,ps_num,DIM)
        prims = Matrix{Cdouble}(undef,ps_num,DIM+2)
        df = Matrix{Cdouble}(undef,vs_nums,NDF)
        level = Vector{Int8}(undef,vs_nums)
        sdata[i] = IBTransferData(midpoints,prims,level,df)
        # sdata[i] = sc_malloc(-1,vnums[i]*(sizeof(Cdouble)*NDF+sizeof(Int8))+ps_num*DIM*sizeof(Cdouble))
        # midpoint_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]),(ps_num,DIM))
        # df_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums[i],NDF))
        # level_buffer = unsafe_wrap(Vector{Int8},Ptr{Int8}(sdata[i]+ps_num*DIM*sizeof(Cdouble)+vnums[i]*sizeof(Cdouble)*NDF),vnums[i])
        # prims_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]+ps_num*DIM*sizeof(Cdouble)+vnums[i]*(sizeof(Cdouble)*NDF+sizeof(Int8))),(ps_num,DIM+2))
        offset = 0
        for j in eachindex(IB_ranks_table[i])
            midpoints[j,:] .= IB_ranks_table[i][j].midpoint
            prims[j,:] .= IB_ranks_table[i][j].prim
        end
        for j in eachindex(IB_ranks_table[i])
            vs_num = IB_ranks_table[i][j].vs_data.vs_num
            df[offset+1:offset+vs_num,:] .= IB_ranks_table[i][j].vs_data.df
            level[offset+1:offset+vs_num] .= IB_ranks_table[i][j].vs_data.level
            offset+=vs_num
        end
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
        ps_num = length(r_vs_nums[i])
        vs_nums = sum(r_vs_nums[i])
        # rdata[i] = sc_malloc(-1,sum(r_vs_nums[i])*(sizeof(Cdouble)*NDF+sizeof(Int8))+(2*DIM+2)*ps_nums*sizeof(Cdouble))
        midpoints = Matrix{Cdouble}(undef,ps_num,DIM)
        prims = Matrix{Cdouble}(undef,ps_num,DIM+2)
        df = Matrix{Cdouble}(undef,vs_nums,NDF)
        level = Vector{Int8}(undef,vs_nums)
        rdata[i] = IBTransferData(midpoints,prims,level,df)
        # rdata[i] = sc_malloc(-1,sum(r_vs_nums[i])*(sizeof(Cdouble)*NDF+sizeof(Int8))+(2*DIM+2)*ps_nums*sizeof(Cdouble))
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
                df,
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
        !isassigned(r_vs_nums,i)&&continue
        ps_offset = 1;vs_offset = 0
        # vnums = sum(r_vs_nums[i])
        # ps_num = length(r_vs_nums[i])
        midpoints = rdata[i].midpoints
        df_buffer = rdata[i].df
        level_buffer = rdata[i].level
        prims = rdata[i].prims
        for j in eachindex(IB_cells)
            for k in eachindex(IB_cells[j].IB_nodes)
                for l in eachindex(IB_cells[j].IB_nodes[k])
                    !isa(IB_cells[j].IB_nodes[k][l],GhostIBNode) && continue
                    vs_num = r_vs_nums[i][ps_offset]
                    node = IB_cells[j].IB_nodes[k][l]
                    node.midpoint = @view(midpoints[ps_offset,:])
                    node.prim = @view(prims[ps_offset,:])
                    node.vs_data.level = @view(level_buffer[vs_offset+1:vs_offset+vs_num])
                    node.vs_data.df = @view(df_buffer[vs_offset+1:vs_offset+vs_num,:])
                    ps_offset+=1;vs_offset+=vs_num
                end
            end
        end
    end
end

function IB_structure_update!(amr::AMR)
    boundary = amr.field.boundary
    IB_vs_nums_exchange!(boundary)
    IB_structure_exchange!(boundary)
    IB_wrap_update!(boundary)
end
function update_quadid!_kernel(ip,data,dp)
    quadid = global_quadid(ip)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    isa(ps_data,InsideSolidData)&&return nothing
    ps_data.quadid = quadid
end
function update_quadid!(info,data)
    AMR_volume_iterate(info, data, P4est_PS_Data, update_quadid!_kernel)
end
function IB_quadid_exchange!(boundary::Boundary)
    rank = MPI.Comm_rank(MPI.COMM_WORLD) 
    solid_cells = boundary.solid_cells
    Numbers = boundary.Numbers
    rbuffer = Vector{Vector{Vector{Cint}}}(undef,length(solid_cells)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Cint}}(undef,length(solid_cells)) # boundaries{points}
    for i in eachindex(solid_cells)
        buffer = Vector{Cint}(undef,length(solid_cells[i].ps_datas))
        for j in eachindex(solid_cells[i].ps_datas)
            buffer[j] = solid_cells[i].ps_datas[j].quadid
        end
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
function IB_solid_cells_assign(amr::AMR{DIM,NDF}) where{DIM,NDF}
    IB_nodes = Vector{Vector{Vector{PS_Data{DIM,NDF}}}}(undef,length(amr.field.boundary.solid_cells))
    for i in eachindex(IB_nodes)
        IB_nodes[i] = Vector{Vector{PS_Data{DIM,NDF}}}(undef,length(amr.field.boundary.solid_cells[i].quadids))
        for j in eachindex(IB_nodes[i])
            IB_nodes[i][j] = PS_Data{DIM,NDF}[]
        end
    end
    trees = amr.field.trees.data
    for i in eachindex(trees)
        for j in eachindex(trees[i])
            ps_data = trees[i][j]
            isa(ps_data,InsideSolidData)&&continue
            if ps_data.bound_enc>0
                solid_cell_indices = solid_cell_index_decoder(ps_data.solid_cell_index)
                for solid_cell_index in solid_cell_indices
                    push!(IB_nodes[ps_data.bound_enc][solid_cell_index],ps_data)
                end
            end
        end
    end
    return IB_nodes
end
function IB_update!(ps4est::P_pxest_t,amr::AMR{DIM,NDF}) where{DIM,NDF} # After partition, IB need to be reinitialized
    IB_nodes = IB_solid_cells_assign(amr)
    # for i in eachindex(IB_nodes)
    #     IB_nodes[i] = Vector{Vector{PS_Data{DIM,NDF}}}(undef,length(solid_cells[i].quadids))
    #     for j in eachindex(IB_nodes[i])
    #         IB_nodes[i][j] = PS_Data{DIM,NDF}[]
    #     end
    # end
    boundary = amr.field.boundary
    solid_cells = boundary.solid_cells
    Numbers,IB_ranks_table,IB_cells = IB_Numbers_ranks(ps4est,IB_nodes,solid_cells)
    # finalize_IB!(boundary.IB_buffer)
    boundary.IB_buffer = init_IBCells(Numbers,solid_cells,IB_ranks_table,IB_cells)
    sort_IB_cells!(amr.global_data,boundary)
    init_solid_cells!(boundary)
    boundary.IB_ranks_table = IB_ranks_table;boundary.IB_cells = IB_cells
end
function colinear_test(p1,p2,p3)
    isapprox(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]),0.;atol=1e-10)
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
function sort_IB_cells!(circle::Circle,config::Configure,solid_cells::SolidCells,::Vector,IB_cells::IBCells)
    for i in eachindex(solid_cells.ps_datas)
        ps_data = solid_cells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        if isa(config.IB_sort,DistanceIBSort)
            basis = [norm(x.midpoint .-aux_point) for x in IB_cells.IB_nodes[i]]
        end
        index = sortperm(basis)
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

function init_solid_cells!(boundary::Boundary{DIM,NDF})where{DIM,NDF}
    solid_cells = boundary.solid_cells
    IB_cells = boundary.IB_cells
    @inbounds for i in eachindex(solid_cells)
        for j in eachindex(solid_cells[i].ps_datas)
            ps_data = solid_cells[i].ps_datas[j]
            vs_data = ps_data.vs_data
            IB_node = first(IB_cells[i].IB_nodes[j])
            IB_vs = IB_node.vs_data
            vs_data.vs_num = IB_vs.vs_num
            vs_data.level = IB_vs.level
            vs_data.weight = IB_vs.weight
			vs_data.midpoint = IB_vs.midpoint
            vs_data.df = IB_vs.df
			vs_data.sdf = zeros(IB_vs.vs_num,NDF,DIM) 
            vs_data.flux = zeros(IB_vs.vs_num,NDF)
        end
    end
end
