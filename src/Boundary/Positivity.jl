function positivity_preserving_ib!(ps_data::PsData{DIM,NDF}, volume, Δt) where {DIM,NDF}
    wb = ps_data.w # Positive part
    if ps_data.bound_enc==0
        return nothing
    else
        vs_data = ps_data.vs_data
        solid_dirs = findall(x->isa(x[1], SolidNeighbor), ps_data.neighbor.data)
        micros = zeros(size(vs_data.df, 1), size(vs_data.df, 2), length(solid_dirs))
        we = zeros(DIM+2)
        for (i, id) in enumerate(solid_dirs)
            solid_neighbor = ps_data.neighbor.data[id][1]
            rot = get_rot(id);
            dir = get_dir(id)
            f_midpoint = copy(ps_data.midpoint)
            f_midpoint[dir] -= 0.5*rot*ps_data.ds[dir]
            nheavi = [x>0.0 for x in rot*@views vs_data.midpoint[:, dir]]
            ndf = @views solid_neighbor.vs_data.flux[nheavi, :]
            there_mid = @views vs_data.midpoint[nheavi, :]
            ndx = [
                f_midpoint[j]-there_mid[i, j]*Δt-solid_neighbor.midpoint[j] for
                i in axes(there_mid, 1), j in axes(there_mid, 2)
            ]
            vn = @views there_mid[:, dir]
            area = rot*reduce(*, ps_data.ds[FAT[DIM-1][dir]])
            micro = @views micros[nheavi, :, i]
            sdf = @views solid_neighbor.vs_data.sdf[nheavi, :, :]
            for j in axes(micro, 2)
                for i in axes(micro, 1)
                    micro[i, j] = @views (ndf[i, j]+dot(ndx[i, :], sdf[i, j, :]))*vn[i]*area
                end
            end
            we .+= @views micro_to_macro(
                micros[:, :, i],
                vs_data.midpoint,
                vs_data.weight,
                vs_data,
            )
        end
        we .*= Δt/volume
        δ = 1e-3
        θρ = we[1] > 0 ? 1.0 : min(1.0, (1-δ)*wb[1]/(abs(we[1])+eps()))
        rub = wb[2:(end-1)];
        rue = we[2:(end-1)]
        eb = wb[end]-(sum(abs2, rub))/(2*wb[1])
        ee = we[end]-dot(rub, rue)/wb[1]
        γ = sum(abs2, rue)/(2*wb[1])
        θe = min(1.0, 2*(1-δ)*eb/(√(ee^2+4*γ*(1-δ)*eb)-ee+eps()))
        θ = min(θρ, θe)
        wb .+= θ*we
        for i in axes(micros, 3)
            vs_data.flux .+= @views θ*micros[:, :, i]
        end
    end
end
