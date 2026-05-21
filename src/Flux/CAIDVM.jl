"""
$(TYPEDSIGNATURES)
"""
function calc_domain_flux(::Type{CAIDVM},here_vs::FaceVsData,face::DomainFace{DIM,NDF,Maxwellian},ka::KA) where{DIM,NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,here_weight,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = ka.kinfo.status.Δt
    vs_data = ps_data.vs_data
    nheavi = [!x for x in heavi]
    there_mid = @views vs_data.midpoint[nheavi,:]
    there_weight = @views vs_data.weight[nheavi]
    there_vn = @views there_mid[:,direction]
    @inbounds @views dx = [midpoint[j]-here_mid[i,j]*Δt-ps_data.midpoint[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
    @inbounds @views df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
    bc = get_bc(face.domain.bc)
    bc[1] = calc_ρw(there_mid,df,bc,here_vn,there_vn,here_weight,there_weight,ka)
    there_df = Matrix{Float64}(undef,size(there_mid,1),NDF)
    for i in axes(there_df,1)
        @views there_df[i,:] .= discrete_maxwell(there_mid[i,:],bc,ka.kinfo)
    end
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [there_df[i,j]*there_vn[i] for i in axes(there_df,1),j in axes(there_df,2)]
    fw = micro_to_macro(here_micro,here_mid,here_weight,vs_data)+micro_to_macro(there_micro,there_mid,there_weight,vs_data)
    return fw,[here_micro,there_micro]
end
"""
$(TYPEDSIGNATURES)
"""
function calc_domain_flux(::Type{CAIDVM},here_vs::FaceVsData,face::DomainFace{DIM,NDF,SuperSonicInflow},ka::KA) where{DIM,NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,here_weight,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = ka.kinfo.status.Δt
    vs_data = ps_data.vs_data
    nheavi = [!x for x in heavi]
    there_mid = @views vs_data.midpoint[nheavi,:]
    there_vn = @views there_mid[:,direction]
    there_weight = @views vs_data.weight[nheavi]
    @inbounds @views dx = [midpoint[j]-here_mid[i,j]*Δt-ps_data.midpoint[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
    @inbounds @views df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
    bc = get_bc(face.domain.bc)
    there_df = Matrix{Float64}(undef,size(there_mid,1),NDF)
    for i in axes(there_df,1)
        @views there_df[i,:] .= discrete_maxwell(there_mid[i,:],bc,ka.kinfo)
    end
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [there_df[i,j]*there_vn[i] for i in axes(there_df,1),j in axes(there_df,2)]
    fw = micro_to_macro(here_micro,here_mid,here_weight,vs_data)+micro_to_macro(there_micro,there_mid,there_weight,vs_data)
    return fw,[here_micro,there_micro]
end
"""
$(TYPEDSIGNATURES)
"""
function calc_domain_flux(::Type{CAIDVM},here_vs::FaceVsData,face::DomainFace{DIM,NDF,UniformOutflow},ka::KA) where{DIM,NDF}
    _,direction,_,_,ps_data = unpack(face)
    heavi,here_weight,here_mid,here_vn,_,_ = unpack(here_vs)
    vs_data = ps_data.vs_data
    nheavi = [!x for x in heavi]
    there_mid = @views vs_data.midpoint[nheavi,:]
    ndf = @views vs_data.df[nheavi,:]
    there_vn = @views there_mid[:,direction]
    there_weight = @views vs_data.weight[nheavi]
    df = @views vs_data.df[heavi,:]
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    fw = micro_to_macro(here_micro,here_mid,here_weight,vs_data)+micro_to_macro(there_micro,there_mid,there_weight,vs_data)
    return fw,[here_micro,there_micro]
end
"""
$(TYPEDSIGNATURES)
"""
function calc_domain_flux(::Type{CAIDVM},here_vs::FaceVsData,face::DomainFace{2,NDF,InterpolatedOutflow},ka::KA) where{NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,here_weight,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = ka.kinfo.status.Δt
    here_ps_mid = ps_data.midpoint
    there_ps_mid = 2.0*midpoint-here_ps_mid
    vs_data = ps_data.vs_data
    nheavi = [!x for x in heavi]
    @inbounds @views begin
        there_mid = vs_data.midpoint[nheavi,:]
        there_sdf = vs_data.sdf[nheavi,:,:]
        there_df = vs_data.df[nheavi,:]+(there_ps_mid[direction]-here_ps_mid[direction])*there_sdf[:,:,direction]
        there_vn = there_mid[:,direction]
        there_weight = vs_data.weight[nheavi]
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx[i,:],there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
        here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
        there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    end
    fw = micro_to_macro(here_micro,here_mid,here_weight,vs_data)+micro_to_macro(there_micro,there_mid,there_weight,vs_data)
    return fw,[here_micro,there_micro]
end
"""
$(TYPEDSIGNATURES)
Compute the flux across a single inner face. For [`CAIDVM`](@ref), both macro and micro fluxes are computed.
"""
function calc_flux(::Type{CAIDVM},here_vs,there_vs,flux_data::Union{FullFace,FluxData},ka::KA{DIM,NDF}) where{DIM,NDF} # without face area and Δt
    _,_,midpoint,here_data,there_data = unpack(flux_data)
    Δt = ka.kinfo.status.Δt
    here_mid = here_vs.midpoint;there_mid = there_vs.midpoint
    here_vn = here_vs.vn;there_vn = there_vs.vn
    here_df = here_vs.df;there_df = there_vs.df
    here_sdf = here_vs.sdf;there_sdf = there_vs.sdf
    here_ps_mid = here_data.midpoint;there_ps_mid = there_data.midpoint
    @inbounds @views begin
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        if there_data.bound_enc<0
            here_micro = [(here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]))*here_vn[i] for i in axes(here_df,1),j in axes(here_df,2)]
            there_micro = [there_df[i,j]*there_vn[i] for i in axes(there_df,1),j in axes(there_df,2)]
        else
            here_micro = positivity_preserving_reconstruct(here_df,here_sdf,here_data.ds,dx,here_vn)
            there_micro = positivity_preserving_reconstruct(there_df,there_sdf,there_data.ds,ndx,there_vn)
        end
    end
    here_weight = here_vs.weight;there_weight = there_vs.weight
    fw = micro_to_macro(here_micro,here_mid,here_weight,here_data.vs_data)+micro_to_macro(there_micro,there_mid,there_weight,here_data.vs_data)
    return fw,[here_micro,there_micro]
end

"""
$(TYPEDSIGNATURES)
Positivity preserving reconstruction (10.1016/j.jcp.2009.12.030).
"""
function positivity_preserving_reconstruct(here_df,here_sdf,here_ds,dx,vn) # positivity preserving reconstruct (10.1016/j.jcp.2009.12.030)
    @views begin
        micro = [(here_df[i,j]+min(abs((here_df[i,j]-eps())/(0.5*dot(here_ds,abs.(here_sdf[i,j,:]))+EPS)),1.)*dot(dx[i,:],here_sdf[i,j,:]))*vn[i] for i in axes(here_df,1),j in axes(here_df,2)]
    end
    return micro
end