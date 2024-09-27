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

function calc_intersect_point(boundaries::Vector{AbstractBoundary},solid_cells::Vector{T}) where{T<:SolidCells}
    aux_points = Vector{Vector{Vector{Float64}}}(undef,length(boundaries))
    for i in eachindex(solid_cells)
        aux_points[i] = calc_intersect_point(boundaries[i],solid_cells[i])
    end
    return aux_points
end
function calc_intersect_point(boundary::AbstractBoundary,solid_cells::SolidCells)
    aux_points = Vector{Vector{Float64}}(undef,length(solid_cells.ps_datas))
    for i in eachindex(solid_cells.ps_datas)
        aux_points[i] = calc_intersect_point(boundary,solid_cells.ps_datas[i].midpoint)
    end
    return aux_points
end
function calc_intersect_point(boundary::Circle,midpoint::AbstractVector)
    @. boundary.radius/norm(midpoint-boundary.center)*(midpoint-boundary.center)+boundary.center
end

function pre_broadcast_boundary_points(solid_cells::Vector{T}) where{T<:SolidCells}
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
function broadcast_boundary_points(solid_cells::Vector{SolidCells{DIM,NDF}},aux_points::Vector{T}) where{DIM,NDF,T<:AuxPoints} # Broadcast solid and aux points
    Numbers = pre_broadcast_boundary_points(solid_cells)
    rbuffer = Vector{Vector{Vector{Float64}}}(undef,length(solid_cells)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Float64}}(undef,length(solid_cells)) # boundaries{points}
    for i in eachindex(solid_cells)
        buffer = Vector{Float64}(undef,DIM*length(solid_cells[i].ps_datas))
        for j in eachindex(solid_cells[i].ps_datas)
            buffer[DIM*(i-1)+1:DIM*i] .= solid_cells[i].ps_datas[j].midpoint
        end
        sbuffer[i] = buffer
    end
    for i in eachindex(solid_cells)
        rbuffer[i] = Vector{Vector{Float64}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
        for j in eachindex(Numbers)
            if j-1==MPI.Comm_rank(MPI.COMM_WORLD) 
                rbuffer[i][j] = sbuffer[i]
            else
                rbuffer[i][j] = Vector{Float64}(undef,DIM*Numbers[j][i])
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
            for k in 1:length(rbuffer[i][j])/DIM
               push!(solid_cells[i].global_midpoints,rbuffer[i][j][DIM*(k-1)+1:DIM*k]) 
            end
        end
    end

    for i in eachindex(solid_cells)
        buffer = Vector{Float64}(undef,DIM*length(solid_cells[i].ps_datas))
        for j in eachindex(solid_cells[i].ps_datas)
            buffer[DIM*(i-1)+1:DIM*i] .= aux_points[i].midpoints[j]
        end
        sbuffer[i] = buffer
    end
    for i in eachindex(solid_cells)
        for j in eachindex(Numbers)
            j-1==MPI.Comm_rank(MPI.COMM_WORLD) && (rbuffer[i][j] = sbuffer[i])
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
            for k in 1:length(rbuffer[i][j])/DIM
               push!(aux_points.midpoints,rbuffer[i][j][DIM*(k-1)+1:DIM*k])
               push!(aux_points.mpi_rank,j-1)
            end
        end
    end
end