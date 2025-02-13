function calc_domain_flux(::DVM,here_vs::Face_VS_Data,face::DomainFace{DIM,NDF,SuperSonicInflow},amr::AMR) where{DIM,NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,_,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = amr.global_data.status.Δt
    vs_data = ps_data.vs_data
    there_mid = @views vs_data.midpoint[.!heavi,:]
    there_vn = @views there_mid[:,direction]
    @inbounds @views dx = [midpoint[j]-here_mid[i,j]*Δt-ps_data.midpoint[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
    @inbounds @views df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
    bc = get_bc(face.domain.bc)
    there_df = Matrix{Float64}(undef,size(there_mid,1),NDF)
    for i in axes(there_df,1)
        @views there_df[i,:] .= discrete_maxwell(there_mid[i,:],bc,amr.global_data)
    end
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [there_df[i,j]*there_vn[i] for i in axes(there_df,1),j in axes(there_df,2)]
    return nothing,[here_micro,there_micro]
end
function calc_domain_flux(::DVM,here_vs::Face_VS_Data,face::DomainFace{DIM,NDF,UniformOutflow},amr::AMR) where{DIM,NDF}
    _,direction,_,_,ps_data = unpack(face)
    heavi,_,_,here_vn,_,_ = unpack(here_vs)
    vs_data = ps_data.vs_data
    there_mid = @views vs_data.midpoint[.!heavi,:]
    ndf = @views vs_data.df[.!heavi,:]
    there_vn = @views there_mid[:,direction]
    df = @views vs_data.df[heavi,:]
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    return nothing,[here_micro,there_micro]
end
function calc_domain_flux(::DVM,here_vs::Face_VS_Data,face::DomainFace{2,NDF,InterpolatedOutflow},amr::AMR) where{NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,_,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = amr.global_data.status.Δt
    here_ps_mid = ps_data.midpoint
    there_ps_mid = 2.0*midpoint-here_ps_mid
    vs_data = ps_data.vs_data
    nheavi = [!x for x in heavi]
    @inbounds @views begin
        there_mid = vs_data.midpoint[nheavi,:]
        there_sdf = vs_data.sdf[nheavi,:,:]
        there_df = vs_data.df[nheavi,:]+(there_ps_mid[direction]-here_ps_mid[direction])*there_sdf[:,:,direction]
        there_vn = there_mid[:,direction]
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx[i,:],there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
        here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
        there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    end
    return nothing,[here_micro,there_micro]
end
function calc_flux(::DVM,here_vs,there_vs,flux_data::Union{FullFace,Flux_Data},amr::AMR{DIM,NDF}) where{DIM,NDF} # without face area and Δt
    _,_,midpoint,here_data,there_data = unpack(flux_data)
    Δt = amr.global_data.status.Δt
    here_mid = here_vs.midpoint;there_mid = there_vs.midpoint
    here_vn = here_vs.vn;there_vn = there_vs.vn
    here_df = here_vs.df;there_df = there_vs.df
    here_sdf = here_vs.sdf;there_sdf = there_vs.sdf
    here_ps_mid = here_data.midpoint;there_ps_mid = there_data.midpoint
    @inbounds @views begin
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx[i,:],there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
        here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
        there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    end
    return nothing,[here_micro,there_micro]
end