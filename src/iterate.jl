function update_gradmax!(amr::AMR)
    global_data = amr.global_data
    global_data.status.gradmax =
        MPI.Allreduce(global_data.status.gradmax, (x,y)->max.(x,y), MPI.COMM_WORLD)
end
function iterate!(amr::AMR)
    time_marching = amr.global_data.config.solver.time_marching
    stable_flag = amr.global_data.status.stable_flag
    (!stable_flag[1]&&stable_flag[2])&&(amr.global_data.status.Δt*=TIME_STEP_CONTRACT_RATIO)
    (stable_flag[1]&&!stable_flag[2])&&(amr.global_data.status.Δt/=TIME_STEP_CONTRACT_RATIO)
    stable_flag[2] = stable_flag[1];stable_flag[1] = true
    iterate!(time_marching,amr)
    stable_check!(amr)
    residual_comm!(amr.global_data)
    amr.global_data.status.ps_adapt_step += 1
    amr.global_data.status.vs_adapt_step += 1
    amr.global_data.status.partition_step += 1
    amr.global_data.status.residual.step += 1
    amr.global_data.status.sim_time+=amr.global_data.status.Δt
end
function stable_check!(amr::AMR)
    MPI.Allreduce(amr.global_data.status.stable_flag[1],&,MPI.COMM_WORLD)
end
function iterate!(::UGKS_Marching,amr::AMR)
    global_data = amr.global_data
    gas = global_data.config.gas
    trees = amr.field.trees
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
            τ = get_τ(prim, gas.μᵣ, gas.ω)
            ps_data.qf .= qf = calc_qf(vs_data, prim_)
            F_ = discrete_maxwell(vs_data.midpoint, prim_, global_data)
            F = discrete_maxwell(vs_data.midpoint, prim, global_data)
            F⁺ = shakhov_part(vs_data.midpoint, F_, prim_, qf, global_data)
            F_ .+= F⁺
            F⁺ = shakhov_part(vs_data.midpoint, F, prim, qf, global_data)
            F .+= F⁺
            f = vs_data.df
            Δt = global_data.status.Δt
            @inbounds @. f = abs((f + 0.5 * Δt * (F / τ + (F_ - f) / τ_)) / (1.0 + 0.5 * Δt / τ) + vs_data.flux / area / (1.0 + 0.5 * Δt / τ))
            vs_data.flux .= 0.0
        end
    end
end
function iterate!(::Euler,amr::AMR{DIM}) where{DIM}
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
            if isa(amr.global_data.config.solver.flux,MicroFlux)
                ps_data.w = calc_w0(ps_data)
            else
                ps_data.w .+= ps_data.flux ./ area
                ps_data.flux .= 0.0
            end
            prim = get_prim(ps_data, global_data)
            # @assert !any(x->isnan(x),prim) `prim NaN!`,ps_data.bound_enc,ps_data.midpoint,prim,ps_data.vs_data.df
            # @assert !any(x->abs(x)>1e10,prim) `prim Too Large!`,ps_data.bound_enc,ps_data.midpoint,prim,ps_data.vs_data.df
            f = vs_data.df
            if 1/prim[end]<1e-3
                @. f = max(f,0.)
                ps_data.w = calc_w0(ps_data)
                prim = get_prim(ps_data,global_data)
                @. f += Δt/area*vs_data.flux
                global_data.status.stable_flag[1]&&(global_data.status.stable_flag[1]=false)
            else
                τ = get_τ(prim, gas.μᵣ, gas.ω)
                ps_data.qf .= qf = calc_qf(vs_data, prim)
                F = discrete_maxwell(vs_data.midpoint, prim, global_data)
                # @assert !any(x->isnan(x),F) `F NaN!`,ps_data.bound_enc,ps_data.midpoint,prim,F
                # @assert !any(x->abs(x)>1e10,F) `F Too Large!`,ps_data.bound_enc,ps_data.midpoint,prim,F
                F .+= shakhov_part(vs_data.midpoint, F, prim, qf, global_data)
                # @assert !any(x->isnan(x),F) `F+ NaN!`,ps_data.bound_enc,ps_data.midpoint,prim,F
                # @assert !any(x->abs(x)>1e10,F) `F+ Too Large!`,ps_data.bound_enc,ps_data.midpoint,prim,F
                # @. f = abs((τ-Δt)/τ*f+Δt/τ*F+Δt/area*vs_data.flux)
                @. f = (τ-Δt)/τ*f+Δt/τ*F+Δt/area*vs_data.flux
            end
            residual_check!(ps_data,prim,global_data)
            ps_data.prim .= prim
            # for fi in f
            #     fi<0. &&(fi=0.)
            # end
            vs_data.flux .= 0.0
        end
    end
end
function residual_check!(ps_data::PS_Data,prim::Vector{Float64},global_data::Global_Data)
    Res = global_data.status.residual
    Res.step%RES_CHECK_INTERVAL!=0&&(return nothing)
    @. Res.sumRes+=(prim-ps_data.prim).^2
    @. Res.sumAvg+=abs(prim)
    return nothing
end
function residual_comm!(global_data::Global_Data)
    Res = global_data.status.residual
    fp = PointerWrapper(global_data.forest.p4est)
    N = fp.global_num_quadrants[]
    Res.step%RES_CHECK_INTERVAL!=0&&(return nothing)
    MPI.Reduce!(Res.sumRes,(x,y)->x.+y,0,MPI.COMM_WORLD)
    MPI.Reduce!(Res.sumAvg,(x,y)->x.+y,0,MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @. Res.residual=sqrt(Res.sumRes*N)/(Res.sumAvg+EPS)
    end
    MPI.Bcast!(Res.residual,0,MPI.COMM_WORLD)
    Res.sumRes.=0.;Res.sumAvg.=0.
end
function check_for_convergence(amr::AMR)
    maximum(amr.global_data.status.residual.residual)<TOLERANCE&&(amr.global_data.status.residual.redundant_step+=1)
    return amr.global_data.status.residual.redundant_step>REDUNDANT_STEPS_NUM
end