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
function broadcast_boundary_midpoints!(boundary_points::Vector{Vector{Vector{Float64}}},image_points::Vector{Vector{Vector{Float64}}},::Global_Data{DIM})where{DIM}
    Numbers = pre_broadcast_boundary_points(boundary_points)
    rbuffer = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Float64}}(undef,length(boundary_points)) # boundaries{points}
    for i in eachindex(boundary_points)
        buffer = Vector{Float64}(undef,2*DIM*length(boundary_points[i]))
        for j in eachindex(boundary_points[i])
            buffer[2*DIM*(i-1)+1:2*DIM*i-DIM] .= boundary_points[i][j]
            buffer[2*DIM*(i-1)+DIM+1:2*DIM*i] .= image_points[i][j]
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
                        ps_data.solid_cell_index = l
                        push!(IB_nodes[k][l],ps_data)
                    end
                end
            end
        end
    end
end

function IB_data_exchange!(amr::AMR)
    boundary = amr.field.boundary
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

function IB_structure_update!(amr::AMR)
    boundary = amr.field.boundary
    IB_vs_nums_exchange!(boundary)
    IB_structure_exchange!(boundary)
    IB_wrap_update!(boundary)
end
function update_quadid!_kernel(ip,data,dp)
    quadid = global_quadid(ip)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
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
    AMR_4est_volume_iterate(ps4est, p_data, update_quadid!)
    IB_quadid_exchange!(boundary::Boundary)
end
function IB_solid_cells_assign(amr::AMR{DIM,NDF}) where{DIM,NDF}
    IB_nodes = Vector{Vector{PS_Data{DIM,NDF}}}(undef,length(amr.field.boundary.solid_cells))
    for i in eachindex(IB_nodes)
        IB_nodes[i] = Vector{PS_Data{DIM,NDF}}(undef,length(amr.field.boundary.solid_cells[i].quadids))
    end
    trees = amr.field.trees.data
    for i in eachindex(trees)
        for j in eachindex(trees[i])
            if trees[i][j].bound_enc>0
                ps_data = trees[i][j]
                push!(IB_nodes[ps_data.bound_enc][ps_data.solid_cell_index],ps_data)
            end
        end
    end
    return IB_nodes
end
function IB_update!(ps4est::P_pxest_t,amr::AMR{DIM,NDF}) where{DIM,NDF} # After partition, IB need to be reinitialized
    IB_nodes = IB_solid_cells_assign(amr)
    for i in eachindex(IB_nodes)
        IB_nodes[i] = Vector{Vector{PS_Data{DIM,NDF}}}(undef,length(solid_cells[i].quadids))
        for j in eachindex(IB_nodes[i])
            IB_nodes[i][j] = PS_Data{DIM,NDF}[]
        end
    end
    boundary = amr.field.boundary
    Numbers,IB_ranks_table,IB_cells = IB_Numbers_ranks(ps4est,IB_nodes,solid_cells)
    finalize_IB!(boundary.IB_buffer)
    boundary.IB_buffer = init_IBCells(Numbers,solid_cells,IB_ranks_table,IB_cells)
    boundary.IB_ranks_table = IB_ranks_table;boundary.IB_cells = IB_cells
end

function sort_IB_cells!(circle::Circle,IB_sort::Symbol,solid_cells::SolidCells,::Vector,IB_cells::IBCells)
    for i in eachindex(solid_cells.ps_datas)
        ps_data = solid_cells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        if IB_sort==:distance
            basis = [norm(x.midpoint .-aux_point) for x in IB_cells.IB_nodes[i]]
        end
        index = sortperm(basis)
        IB_cells.IB_nodes[i] = IB_cells.IB_nodes[i][index]
    end
end
function sort_IB_cells!(global_data::Global_Data,boundary::Boundary)
    IBs = global_data.config.IB
    solid_cells = boundary.solid_cells
    aux_points = boundary.aux_points
    IB_cells = boundary.IB_cells
    for i in eachindex(IBs)
        sort_IB_cells!(IBs[i],global_data.config.IB_sort,solid_cells[i],aux_points[i],IB_cells[i])
    end
end