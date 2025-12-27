function update_gradmax!(amr::AMR{DIM}) where{DIM}
    global_data = amr.global_data
    update_sw!(amr)
    global_data.status.gradmax =
        MPI.Allreduce(global_data.status.gradmax, (x,y)->max.(x,y), MPI.COMM_WORLD)
    !global_data.config.solver.AMR_PS_SMOOTH&&return nothing
    ps_datas = amr.ghost.ghost_wrap
    for ps_data in ps_datas
        isa(ps_data,AbstractInsideSolidData)&&continue
        vs_data = ps_data.vs_data
        for i in 1:DIM
            ps_data.sw[:,i].=@views calc_w0(vs_data.midpoint,vs_data.sdf[:,:,i],vs_data.weight,amr.global_data)
        end
    end
    for tree in amr.field.trees.data
        for ps_data in tree
            isa(ps_data, InsideSolidData)&&continue
            ps_data.bound_enc!=0&&continue
            ps_data.flux.=0.
            for neighbors in ps_data.neighbor.data
                for neighbor in neighbors
                    isnothing(neighbor)&&continue
                    for j in axes(ps_data.sw,1)
                        for k in axes(ps_data.sw,2)
                            ps_data.flux[j] = max(abs(ps_data.sw[j,k]),max(abs(neighbor.sw[j,k]),abs(ps_data.flux[j]))) # Temporarily used to store local_grad_max
                        end
                    end
                end
            end
        end
    end
    return nothing
end
function iterate!(amr::AMR;buffer_steps = 0,i = typemax(Int))
    time_marching = amr.global_data.config.solver.time_marching
    iterate!(time_marching,amr)
    residual_comm!(amr.global_data)
    amr.global_data.status.ps_adapt_step += 1
    amr.global_data.status.vs_adapt_step += 1
    amr.global_data.status.partition_step += 1
    amr.global_data.status.residual.step += 1
    amr.global_data.status.sim_time+=amr.global_data.status.Δt
end
function stable_check!(amr::AMR) # deprecated 
    MPI.Allreduce(amr.global_data.status.stable_flag[1],&,MPI.COMM_WORLD)
end
function iterate!(::UGKS_Marching,amr::AMR)
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
                global_data.status.Δt=min(TIME_STEP_CONTRACT_RATIO,global_data.config.gas.Kn)*global_data.status.Δt_ξ
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
function iterate!(::CAIDVM_Marching,amr::AMR;buffer_steps = 0, i = typemax(Int))
    global_data = amr.global_data
    gas = global_data.config.gas
    trees = amr.field.trees
    Δt = global_data.status.Δt
    global_data.status.Δt = i>buffer_steps ? global_data.status.Δt_ξ : global_data.status.Δt_ξ/buffer_steps
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
            # if 1/prim_c[end]<1e-6
            #     f.+= Δt/area*vs_data.flux # Convection first
            #     @. f = max(f,0.)
            #     ps_data.w = calc_w0(ps_data)
            #     prim = prim_c = get_prim(ps_data,global_data)
            #     τ = global_data.config.gas.Kn
            #     # global_data.status.Δt=min(TIME_STEP_CONTRACT_RATIO,global_data.config.gas.Kn)*global_data.status.Δt_ξ
            # else
                f.+= Δt/area*vs_data.flux # Convection first
                w = calc_w0(vs_data.midpoint,f,vs_data.weight,global_data)
                prim = get_prim(w,global_data)
                # if 1/prim[end]>5||1/prim[end]<0.
                #     if ps_data.bound_enc>0
                #         sdirs = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
                #         down_id = []
                #         for id in sdirs
                #             sn = ps_data.neighbor.data[id][1]
                #             normal = sn.normal
                #             dir1 = get_dir(id);rot1 = get_rot(id)
                #             append!(down_id,findall(x->dot(x,normal)<0&&x[dir1]*rot1>0,eachrow(vs_data.midpoint)))
                #         end
                #     end
                #     df_i,minid = findmin(vs_data.df[:,1]);df_di,_ = findmin(@views vs_data.df[down_id,1]);
                #     df_j,maxid = findmax(vs_data.df[:,1]);df_dj,_ = findmax(@views vs_data.df[down_id,1]);
                #     @show vs_data.df down_id df_i df_di df_j df_dj vs_data.midpoint[minid,:] vs_data.midpoint[maxid,:] ps_data.midpoint
                #     throw(`instale error!`)
                # end
                # if 1/prim[end]<1e-6
                #     @. f = max(f,0.)
                #     ps_data.w = calc_w0(ps_data)
                #     prim = prim_c = get_prim(ps_data,global_data)
                #     τ = global_data.config.gas.Kn
                #     # global_data.status.Δt=min(TIME_STEP_CONTRACT_RATIO,global_data.config.gas.Kn)*global_data.status.Δt_ξ
                # else
                    τ = get_τ(prim_c, gas.μᵣ, gas.ω) # τ^{n+1}
                # end
            # end
            # if any(x->(abs(x)>1e10),prim_c)
            #     types = [typeof(x[1]) for x in ps_data.neighbor.data]
            #     @show types f ps_data.w prim_c ps_data.bound_enc
            #     throw(`iterate singular!`)
            # end
            # if prim_c[end]<0
            #     types = [typeof(x[1]) for x in ps_data.neighbor.data]
            #     midpoints = [ps_data.neighbor.data[i][1].midpoint[j] for i in 1:6,j in 1:3]
            #     @show ps_data.bound_enc prim_c prim types ps_data.midpoint midpoints
            # end
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
            # global_data.status.Δt=min(global_data.status.Δt,macro_cons_time_step(ps_data,global_data))
        end
    end
    # Δt_comm!(global_data)
end
function iterate!(::Euler,amr::AMR{DIM}) where{DIM}
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
function Δt_comm!(global_data::Global_Data)
    global_data.status.Δt = MPI.Allreduce(global_data.status.Δt, (x,y)->min(x,y), MPI.COMM_WORLD)
end
function macro_cons_time_step(ps_data,global_data::Global_Data{DIM}) where{DIM}
    dt = Inf;CFL = global_data.config.solver.CFL
    prim = ps_data.prim
    for i in 2:DIM+1
        c = √(prim[i]^2+0.5/prim[end])
        dt = min(ps_data.ds[i-1]/(c/(erf(c*√prim[end])+EPS)),dt)
    end
    return CFL*dt
end