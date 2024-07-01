function update_Δt!(DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    trees = DVM_data.trees
    @inbounds @simd for i in eachindex(trees.data)
        @inbounds @simd for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            U = global_data.quadrature[2] - 10 / 56
            V = global_data.quadrature[4] - 10 / 56
            Δt = 1.0
            Δt = min(Δt, global_data.gas.CFL * min(ps_data.ds[1] / U, ps_data.ds[2] / V))
            global_data.gas.Δt = Δt
        end
    end
    min_Δt = MPI.Allreduce(global_data.gas.Δt, min, MPI.COMM_WORLD)
    global_data.gas.Δt = min_Δt
    global_data.gas.sim_time += global_data.gas.Δt
end
function update_gradmax!(DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    global_data.gradmax = MPI.Allreduce(global_data.gradmax, max, MPI.COMM_WORLD)
end
function update_volume!(DVM_data::DVM_Data)
    gas = DVM_data.global_data.gas
    trees = DVM_data.trees
    @inbounds @simd for i in eachindex(trees.data)
        @inbounds @simd for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            vs_data = ps_data.vs_data
            h = @views vs_data.df[:, 1]
            b = @views vs_data.df[:, 2]
            u = @views vs_data.midpoint[:, 1]
            v = @views vs_data.midpoint[:, 2]
            prim_ = get_prim(ps_data.w, gas.γ)
            τ_ = get_τ(prim_, gas.μᵣ, gas.ω)
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux ./ area
            # if ps_data.midpoint == [0.90625, 0.90625]
            #     @show ps_data.flux ./area
            # end
            ps_data.flux .= 0.0
            ps_data.prim .= prim = get_prim(ps_data.w, gas.γ)
            τ = get_τ(prim, gas.μᵣ, gas.ω)
            ps_data.qf .= qf = calc_qf(vs_data.df, prim_, vs_data.midpoint, vs_data.weight)
            H_, B_ = discrete_maxwell(u, v, prim_, gas.K)
            H, B = discrete_maxwell(u, v, prim, gas.K)
            H⁺, B⁺ = shakhov_part(H_, B_, prim_, u, v, qf, gas.Pr, gas.K)
            H_ .+= H⁺
            B_ .+= B⁺
            H⁺, B⁺ = shakhov_part(H, B, prim, u, v, qf, gas.Pr, gas.K)
            H .+= H⁺
            B .+= B⁺
            vs_data.df[:, 1] .=
                @. (h + 0.5 * gas.Δt * (H / τ + (H_ - h) / τ_)) / (1.0 + 0.5 * gas.Δt / τ)
            vs_data.df[:, 2] .=
                @. (b + 0.5 * gas.Δt * (B / τ + (B_ - b) / τ_)) / (1.0 + 0.5 * gas.Δt / τ)
            @. vs_data.df += vs_data.flux / area / (1.0 + 0.5 * gas.Δt / τ)
            vs_data.flux .= 0.0
            # @show sum(vs_data.df[:,1].*vs_data.weight)
        end
    end
end
