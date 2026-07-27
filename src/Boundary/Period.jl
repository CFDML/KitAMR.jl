function periodic_ghost_cell(
    midpoint::Vector,
    ps_data::AbstractPsData{DIM,NDF},
) where {DIM,NDF}
    # PsData{DIM,NDF}(ps_data.quadid,ps_data.bound_enc,ps_data.solid_cell_index,
    #     ps_data.ds,midpoint,ps_data.qf,ps_data.w,ps_data.sw,ps_data.prim,ps_data.flux,ps_data.vs_data,
    #     ps_data.neighbor
    # )
    PsData(ps_data; midpoint)
end
function periodic_ghost_cell(
    midpoint::Vector,
    ps_data::AbstractGhostPsData{DIM,NDF},
) where {DIM,NDF}
    GhostPsData{DIM,NDF}(
        ps_data.owner_rank,
        ps_data.quadid,
        ps_data.bound_enc,
        ps_data.ds,
        midpoint,
        ps_data.w,
        ps_data.sw,
        ps_data.vs_data,
    )
end
function periodic_ghost_cell(
    midpoints::Vector{Vector{Float64}},
    ps_datas::Vector{NeighborQuad{DIM,NDF}},
) where {DIM,NDF}
    datas = Vector{AbstractPsData{DIM,NDF}}(undef, length(midpoints))
    for i in eachindex(midpoints)
        datas[i] = periodic_ghost_cell(midpoints[i], ps_datas[i])
    end
    return datas
end
