function update_Δt!(amr::AMR{2,NDF}) where {NDF}
    global_data = amr.global_data
    trees = amr.field.trees
    quadrature = global_data.config.quadrature
    U =
        quadrature[2] -
        (quadrature[2] - quadrature[1]) / global_data.config.vs_trees_num[1]/2^global_data.config.solver.AMR_VS_MAXLEVEL / 2
    V =
        quadrature[4] -
        (quadrature[4] - quadrature[3]) / global_data.config.vs_trees_num[2]/2^global_data.config.solver.AMR_VS_MAXLEVEL / 2
    @inbounds @simd for i in eachindex(trees.data)
        @inbounds @simd for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            Δt = 1.0
            Δt = min(Δt, global_data.config.solver.CFL * min(ps_data.ds[1] / U, ps_data.ds[2] / V))
            global_data.status.Δt = Δt
        end
    end
    min_Δt = MPI.Allreduce(global_data.status.Δt, min, MPI.COMM_WORLD)
    global_data.status.Δt = min_Δt
    global_data.status.sim_time += global_data.status.Δt
end
function update_Δt!(amr::AMR{3,NDF}) where {NDF}
    global_data = amr.global_data
    trees = amr.field.trees
    quadrature = global_data.config.quadrature
    U =
        quadrature[2] -
        (quadrature[2] - quadrature[1]) / global_data.config.vs_trees_num[1]/2^global_data.config.solver.AMR_VS_MAXLEVEL / 2
    V =
        quadrature[4] -
        (quadrature[4] - quadrature[3]) / global_data.config.vs_trees_num[2]/2^global_data.config.solver.AMR_VS_MAXLEVEL/ 2
    W =
        quadrature[6] -
        (quadrature[6] - quadrature[5]) / global_data.config.vs_trees_num[3]/2^global_data.config.solver.AMR_VS_MAXLEVEL / 2
    @inbounds @simd for i in eachindex(trees.data)
        @inbounds @simd for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            Δt = 1.0
            Δt = min(
                Δt,
                global_data.config.solver.CFL *
                min(min(ps_data.ds[1] / U, ps_data.ds[2] / V), ps_data.ds[3] / W),
            )
            global_data.status.Δt = Δt
        end
    end
    min_Δt = MPI.Allreduce(global_data.status.Δt, min, MPI.COMM_WORLD)
    global_data.status.Δt = min_Δt
    global_data.status.sim_time += global_data.status.Δt
end
function update_gradmax!(amr::AMR)
    global_data = amr.global_data
    global_data.status.gradmax =
        MPI.Allreduce(global_data.status.gradmax, max, MPI.COMM_WORLD)
end
function update_df!(
    vs_data::VS_Data,
    F::AbstractArray,
    F_::AbstractArray,
    τ::Real,
    τ_::Real,
    area::Real,
    global_data::Global_Data,
)
    f = vs_data.df
    Δt = global_data.status.Δt
    @inbounds @. f = (f + 0.5 * Δt * (F / τ + (F_ - f) / τ_)) / (1.0 + 0.5 * Δt / τ)
    @inbounds @. f += vs_data.flux / area / (1.0 + 0.5 * Δt / τ)
end
function update_volume!(amr::AMR)
    global_data = amr.global_data
    gas = global_data.config.gas
    trees = amr.field.trees
    @inbounds @simd for i in eachindex(trees.data)
        @inbounds @simd for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            vs_data = ps_data.vs_data
            prim_ = get_prim(ps_data, global_data)
            τ_ = get_τ(prim_, gas.μᵣ, gas.ω)
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux ./ area
            ps_data.flux .= 0.0
            ps_data.prim .= prim = get_prim(ps_data, global_data)
            τ = get_τ(prim, gas.μᵣ, gas.ω)
            ps_data.qf .= qf = calc_qf(vs_data, prim_)
            F_ = discrete_maxwell(vs_data.midpoint, prim_, global_data)
            F = discrete_maxwell(vs_data.midpoint, prim, global_data)
            F⁺ = shakhov_part(vs_data.midpoint, F_, prim_, qf, global_data)
            F_ .+= F⁺
            F⁺ = shakhov_part(vs_data.midpoint, F, prim, qf, global_data)
            F .+= F⁺
            update_df!(vs_data, F, F_, τ, τ_, area, global_data)
            vs_data.flux .= 0.0
        end
    end
end
