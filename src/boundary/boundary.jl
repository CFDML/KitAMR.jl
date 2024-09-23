function boundary_interpolation_points(amr::AMR,nodes::AbstractMatrix,r::Real) # Search fluid points inside the detect radius r
    interp_points = Vector{Float64}
    datas = amr.field.trees.data
    for i in eachindex(data)
        for j in eachindex(data[i])
            ps_data = data[i][j]

        end
    end
end

function boundary_ghost_cell() # Search ghost cells, that is those inside solid domain and immediately adjacent the boundary.
end

function boundary_communicate()
end
function init_IB!(global_data::Global_Data{DIM}) where{DIM}
    length(global_data.config.boundary)==2*DIM && return nothing
    boundary = global_data.config.boundary
    for i in eachindex(boundary)
        init_solid_cells!(boundary[i],)
    end
end