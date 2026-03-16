"""
$(TYPEDSIGNATURES)
Outer function of collision and time marching.
"""
function iterate!(amr::KitAMR_Data)
    time_marching = amr.global_data.config.solver.time_marching
    iterate!(time_marching,amr)
    residual_comm!(amr.global_data)
    amr.global_data.status.ps_adapt_step += 1
    amr.global_data.status.vs_adapt_step += 1
    amr.global_data.status.partition_step += 1
    amr.global_data.status.residual.step += 1
    amr.global_data.status.sim_time+=amr.global_data.status.Δt
    return nothing
end
function iterate!(::Type{UGKS_Marching},amr::KitAMR_Data)
    global_data = amr.global_data
    gas = global_data.config.gas
    trees = amr.field.trees
    Δt = global_data.status.Δt
    global_data.status.Δt = global_data.status.Δt_ξ
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            prim_ = get_prim(ps_data, global_data)
            τ_ = get_τ(prim_, gas.μᵣ, gas.ω)
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux ./ area
            ps_data.flux .= 0.0
            ps_data.prim .= prim = get_prim(ps_data, global_data)
            f = vs_data.df
            if 1/prim[end]<1e-6
                @. f = max(f,0.)
                ps_data.w = calc_w0(ps_data)
                prim = get_prim(ps_data,global_data)
                τ = global_data.config.gas.Kn
            else
                τ = get_τ(prim, gas.μᵣ, gas.ω)
            end
            residual_check!(ps_data,prim_,global_data)
            ps_data.qf .= qf = calc_qf(vs_data, prim_)
            F_ = discrete_maxwell(vs_data.midpoint, prim_, global_data)
            F = discrete_maxwell(vs_data.midpoint, prim, global_data)
            F⁺ = shakhov_part(vs_data.midpoint, F_, prim_, qf, global_data)
            F_ .+= F⁺
            F⁺ = shakhov_part(vs_data.midpoint, F, prim, qf, global_data)
            F .+= F⁺
            f = vs_data.df
            @inbounds @. f = abs((f + 0.5 * Δt * (F / τ + (F_ - f) / τ_)) / (1.0 + 0.5 * Δt / τ) + vs_data.flux / area / (1.0 + 0.5 * Δt / τ))
            vs_data.flux .= 0.0
        end
    end
    Δt_comm!(global_data)
end
# Conserved Adaptive Implicit DVM (CAIDVM)
"""
$(TYPEDSIGNATURES)
Iteration for Conserved Adaptive Implicit DVM (CAIDVM).
"""
function iterate!(::Type{CAIDVM_Marching},amr::KitAMR_Data)
    global_data = amr.global_data
    gas = global_data.config.gas
    trees = amr.field.trees
    Δt = global_data.status.Δt
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux .*Δt / area # Macroscopic update
            prim_c = get_prim(ps_data, global_data) # Conserved macroscopic variables
            f = vs_data.df
            f.+= Δt/area*vs_data.flux # Convection first
            w = calc_w0(vs_data.midpoint,f,vs_data.weight,global_data)
            prim = get_prim(w,global_data)
            τ = get_τ(prim_c, gas.μᵣ, gas.ω) # τ^{n+1}
            F_c = discrete_maxwell(vs_data.midpoint, prim_c, global_data)
            F = discrete_maxwell(vs_data.midpoint, prim, global_data)
            @. f += F_c-F # Conservation correction
            ps_data.qf .= qf = calc_qf(vs_data, prim_c) # Heatflux after convection
            F_c .+= shakhov_part(vs_data.midpoint, F_c, prim_c, qf, global_data) # g^{S,n+1}
            # Collision process
            f .*= τ/(τ+Δt)
            @. f += Δt/(τ+Δt)*F_c
            residual_check!(ps_data,prim_c,global_data)
            ps_data.prim .= prim_c
            ps_data.flux .= 0.0
            vs_data.flux .= 0.0
        end
    end
    return nothing
end
function iterate!(::Type{Euler},amr::KitAMR_Data{DIM}) where{DIM}
    global_data = amr.global_data
    gas = global_data.config.gas
    trees = amr.field.trees
    Δt = global_data.status.Δt
    global_data.status.Δt = global_data.status.Δt_ξ
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            area = reduce(*, ps_data.ds)
            if isa(amr.global_data.config.solver.flux,MicroFlux)
                ps_data.w = calc_w0(ps_data)
            else
                ps_data.w .+= ps_data.flux ./ area
                ps_data.flux .= 0.0
            end
            prim = get_prim(ps_data, global_data)
            f = vs_data.df
            τ = get_τ(prim, gas.μᵣ, gas.ω)
            ps_data.qf .= qf = calc_qf(vs_data, prim)
            F = discrete_maxwell(vs_data.midpoint, prim, global_data)
            F .+= shakhov_part(vs_data.midpoint, F, prim, qf, global_data)
            @. f = (τ-Δt)/τ*f+Δt/τ*F+Δt/area*vs_data.flux
            residual_check!(ps_data,prim,global_data)
            ps_data.prim .= prim
            vs_data.flux .= 0.0
        end
    end
end