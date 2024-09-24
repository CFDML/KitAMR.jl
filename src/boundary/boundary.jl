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