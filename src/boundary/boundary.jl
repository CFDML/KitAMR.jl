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
    aux_points = Vector{Vector{Float64}}(undef,length(solid_cells.ps_datas))
    for i in eachindex(solid_midpoints)
        aux_points[i] = calc_intersect_point(boundary,solid_midpoints[i])
    end
    return aux_points
end
function calc_intersect_point(boundary::Circle,midpoint::AbstractVector)
    @. boundary.radius/norm(midpoint-boundary.center)*(midpoint-boundary.center)+boundary.center
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
end
function broadcast_boundary_midpoints!(boundary_points::Vector{Vector{Vector{Float64}}})
    Numbers = pre_broadcast_aux_points(boundary_points)
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
        for j in eachindex(Numbers)
            if j-1 == MPI.Comm_rank(MPI.COMM_WORLD)
                append!(boundary_points_global[i],boundary_points[i])
            else
                for k in 1:length(rbuffer[i][j])/DIM
                    push!(boundary_points_global[i],rbuffer[i][j][DIM*(k-1)+1:DIM*k]) 
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
            push!(solid_cells[i].quadids,rbuffer[i][j])
        end
        sort!(solid_cells[i].quadids)
    end
end

function IB_flag(aux_point::AbstractVector,midpoint::AbstractVector,ds::AbstractVector,r::Real)
    DIM = length(aux_point)
    for i in eachindex(2^DIM)
        norm(midpoint.+0.5*ds.*RMT[DIM][i].-aux_point)<=r && return true # any corner cross boundary?
    end
    return false
end
function IB_flag(boundaries::Vector{AbstractBoundary},aux_points::Vector{Vector{Float64}},midpoint::AbstractVector,ds::AbstractVector)
    for i in eachindex(boundaries)
        r = boundaries[i].search_radius
        for j in eachindex(aux_points[i])
            IB_flag(aux_points[i][j],midpoint,ds,r) && return true
        end
    end
    return false
end

function search_IB!(IB_nodes::Vector{Vector{Vector{PS_Data{DIM,NDF}}}},aux_points::Vector,trees::PS_Trees{DIM,NDF},global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    ps_datas = trees.data
    boundaries = global_data.config.IB
    for i in eachindex(ps_datas)
        for j in eachindex(ps_datas[i])
            isa(ps_data,InsideSolidData) && continue
            ps_data = ps_datas[i][j]
            for k in eachindex(boundaries)
                r = boundaries[k].search_radius
                for l in eachindex(aux_points[k])
                    IB_flag(aux_points[k][l],ps_data.midpoint,ps_data.ds,r)&&push!(IB_nodes[k][l],ps_data)
                end
            end
        end
    end
end