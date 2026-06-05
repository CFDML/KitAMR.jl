"""
$(TYPEDSIGNATURES)
Outer function of collision and time marching.
"""
function iterate!(p4est::P_pxest_t,ka::KA)
    time_marching = ka.kinfo.config.solver.time_marching
    iterate!(time_marching,ka)
    residual_comm!(ka.kinfo)
    data_exchange!(p4est, ka)
    ka.kinfo.status.step += 1
    ka.kinfo.status.ps_adapt_step += 1
    ka.kinfo.status.vs_adapt_step += 1
    ka.kinfo.status.partition_step += 1
    ka.kinfo.status.sim_time+=ka.kinfo.status.Δt
    return nothing
end
"""
$(TYPEDSIGNATURES)
Shrink the time step `status.Δt` (starting from the grid/CFL step `Δt_ξ`) so the simulation
lands exactly on the next output target: the next animation frame time (integer multiples of
`output.anim_dt`, when animation is enabled) and the termination time `solver.max_sim_time`.
Call once per step before [`flux!`](@ref). With the default `max_sim_time=Inf` and animation
off, this is a no-op (`Δt = Δt_ξ`). Note: only the time-accurate marching schemes
(`CIP_Marching`, `CAIDVM_Marching`) consume `status.Δt` directly; `UGKS`/`Euler` reset it to
`Δt_ξ` mid-step and are meant for steady runs.
"""
function limit_Δt!(ka::KA)
    status = ka.kinfo.status
    output = ka.kinfo.config.output
    Δt = status.Δt_ξ
    if output.anim_dt > 0.
        Δt_anim = (output.anim_index+1)*output.anim_dt - status.sim_time
        Δt_anim > 0. && (Δt = min(Δt, Δt_anim))
    end
    Δt_end = ka.kinfo.config.solver.max_sim_time - status.sim_time
    Δt_end > 0. && (Δt = min(Δt, Δt_end))
    status.Δt = Δt
    return nothing
end
"""
$(TYPEDSIGNATURES)
Return `true` once `sim_time` has reached `solver.max_sim_time` (within a small relative
tolerance). Use as the main-loop termination condition. Always `false` for the default
`max_sim_time=Inf` (i.e. terminate only by convergence).
"""
function reached_max_time(ka::KA)
    return ka.kinfo.status.sim_time >= ka.kinfo.config.solver.max_sim_time*(1-1e-12)
end
function iterate!(::Type{UGKS_Marching},ka::KA)
    kinfo = ka.kinfo
    gas = kinfo.config.gas
    trees = ka.kdata.field.trees
    Δt = kinfo.status.Δt
    kinfo.status.Δt = kinfo.status.Δt_ξ
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            prim_ = get_prim(ps_data, kinfo)
            τ_ = get_τ(prim_, gas.μᵣ, gas.ω)
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux ./ area
            ps_data.flux .= 0.0
            ps_data.prim .= prim = get_prim(ps_data, kinfo)
            f = vs_data.df
            if 1/prim[end]<1e-6
                @. f = max(f,0.)
                ps_data.w = calc_w0(ps_data)
                prim = get_prim(ps_data,kinfo)
                τ = kinfo.config.gas.Kn
            else
                τ = get_τ(prim, gas.μᵣ, gas.ω)
            end
            residual_check!(ps_data,prim_,kinfo)
            ps_data.qf .= qf = calc_qf(vs_data, prim_)
            F_ = discrete_maxwell(vs_data.midpoint, prim_, kinfo)
            F = discrete_maxwell(vs_data.midpoint, prim, kinfo)
            F⁺ = shakhov_part(vs_data.midpoint, F_, prim_, qf, kinfo)
            F_ .+= F⁺
            F⁺ = shakhov_part(vs_data.midpoint, F, prim, qf, kinfo)
            F .+= F⁺
            f = vs_data.df
            @inbounds @. f = abs((f + 0.5 * Δt * (F / τ + (F_ - f) / τ_)) / (1.0 + 0.5 * Δt / τ) + vs_data.flux / area / (1.0 + 0.5 * Δt / τ))
            vs_data.flux .= 0.0
        end
    end
    Δt_comm!(kinfo)
end
# Conserved Adaptive Implicit DVM (CAIDVM)
"""
$(TYPEDSIGNATURES)
Iteration for Conserved Adaptive Implicit DVM (CAIDVM).
"""
function iterate!(::Type{CAIDVM_Marching},ka::KA)
    kinfo = ka.kinfo
    gas = kinfo.config.gas
    trees = ka.kdata.field.trees
    Δt = kinfo.status.Δt
    for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux .*Δt / area # Macroscopic update
            prim_c = get_prim(ps_data, kinfo) # Conserved macroscopic variables
            f = vs_data.df
            f.+= Δt/area*vs_data.flux # Convection first
            w = calc_w0(vs_data.midpoint,f,vs_data.weight,kinfo)
            prim = get_prim(w,kinfo)
            τ = get_τ(prim_c, gas.μᵣ, gas.ω) # τ^{n+1}
            F_c = discrete_maxwell(vs_data.midpoint, prim_c, kinfo)
            F = discrete_maxwell(vs_data.midpoint, prim, kinfo)
            @. f += F_c-F # Conservation correction
            ps_data.qf .= qf = calc_qf(vs_data, prim_c) # Heatflux after convection
            F_c .+= shakhov_part(vs_data.midpoint, F_c, prim_c, qf, kinfo) # g^{S,n+1}
            # Collision process
            f .*= τ/(τ+Δt)
            @. f += Δt/(τ+Δt)*F_c
            residual_check!(ps_data,prim_c,kinfo)
            ps_data.prim .= prim_c
            ps_data.flux .= 0.0
            vs_data.flux .= 0.0
        end
    end
    return nothing
end
function iterate!(::Type{Euler},ka::KA{DIM}) where{DIM}
    kinfo = ka.kinfo
    gas = kinfo.config.gas
    trees = ka.kdata.field.trees
    Δt = kinfo.status.Δt
    kinfo.status.Δt = kinfo.status.Δt_ξ
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            area = reduce(*, ps_data.ds)
            if isa(ka.kinfo.config.solver.flux,MicroFlux)
                ps_data.w = calc_w0(ps_data)
            else
                ps_data.w .+= ps_data.flux .*Δt / area
                ps_data.flux .= 0.0
            end
            prim = get_prim(ps_data, kinfo)
            f = vs_data.df
            τ = get_τ(prim, gas.μᵣ, gas.ω)
            ps_data.qf .= qf = calc_qf(vs_data, prim)
            F = discrete_maxwell(vs_data.midpoint, prim, kinfo)
            F .+= shakhov_part(vs_data.midpoint, F, prim, qf, kinfo)
            f .= (f+Δt/area*vs_data.flux)*τ/(τ+Δt)+Δt/(τ+Δt)*F
            residual_check!(ps_data,prim,kinfo)
            ps_data.prim .= prim
            vs_data.flux .= 0.0
        end
    end
end