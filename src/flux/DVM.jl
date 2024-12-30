function calc_domain_flux(::DVM,here_vs::Face_VS_Data,face::DomainFace{2,NDF,SuperSonicInflow},amr::AMR) where{NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,_,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = amr.global_data.status.Δt
    vs_data = ps_data.vs_data
    there_mid = @views vs_data.midpoint[.!heavi,:]
    there_vn = @views there_mid[:,direction]
    # there_weight = @views vs_data.weight[.!heavi]
    @inbounds @views dx = [midpoint[j]-here_mid[i,j]*Δt-ps_data.midpoint[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
    @inbounds @views df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
    bc = get_bc(face.domain.bc)
    there_df = Matrix{Float64}(undef,size(there_mid,1),NDF)
    for i in axes(there_df,1)
        @views there_df[i,:] .= discrete_maxwell(there_mid[i,:],bc,amr.global_data)
    end
    # flux = Vector{Float64}(undef,4)
    # @views @inbounds flux[1] = sum(here_weight.*here_vn.*here_df[:,1])+sum(there_weight.*there_vn.*there_df[:,1])
    # @views @inbounds flux[2] = sum(here_weight.*here_vn.*df[:,1].*here_mid[:,1])+sum(there_weight.*there_vn.*there_df[:,1].*there_mid[:,1])
    # @views @inbounds flux[3] = sum(here_weight.*here_vn.*df[:,1].*here_mid[:,2])+sum(there_weight.*there_vn.*there_df[:,1].*there_mid[:,2])
    # @views @inbounds flux[end] = 0.5*(sum(here_weight.*here_vn.*(sum(x->x^2,here_mid,dims = 2).*df[:,1]+df[:,2]))+
    #     sum(there_weight.*there_vn.*(sum(x->x^2,there_mid,dims = 2).*there_df[:,1]+there_df[:,2])))
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [there_df[i,j]*there_vn[i] for i in axes(there_df,1),j in axes(there_df,2)]
    return nothing,[Δt*here_micro,Δt*there_micro]
end
function calc_domain_flux(::DVM,here_vs::Face_VS_Data,face::DomainFace{2,NDF,UniformOutflow},amr::AMR) where{NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,_,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = amr.global_data.status.Δt
    here_ps_mid = ps_data.midpoint
    vs_data = ps_data.vs_data
    there_mid = @views vs_data.midpoint[.!heavi,:]
    ndf = @views vs_data.df[.!heavi,:]
    there_vn = @views there_mid[:,direction]
    # there_weight = @views vs_data.weight[.!heavi]
    @inbounds @views dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
    @inbounds @views df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
    # dx_there = midpoint-here_ps_mid
    # flux = Vector{Float64}(undef,4)
    # @views @inbounds flux[1] = sum(here_weight.*here_vn.*here_df[:,1])+sum(there_weight.*there_vn.*ndf[:,1])
    # @views @inbounds flux[2] = sum(here_weight.*here_vn.*df[:,1].*here_mid[:,1])+sum(there_weight.*there_vn.*ndf[:,1].*there_mid[:,1])
    # @views @inbounds flux[3] = sum(here_weight.*here_vn.*df[:,1].*here_mid[:,2])+sum(there_weight.*there_vn.*ndf[:,1].*there_mid[:,2])
    # @views @inbounds flux[end] = 0.5*(sum(here_weight.*here_vn.*(sum(x->x^2,here_mid,dims = 2).*df[:,1]+df[:,2]))+
    #     sum(there_weight.*there_vn.*(sum(x->x^2,there_mid,dims = 2).*ndf[:,1]+ndf[:,2])))
    @inbounds here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
    @inbounds there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    return nothing,[Δt*here_micro,Δt*there_micro]
end
function calc_flux(::DVM,here_vs,there_vs,flux_data::Union{FullFace,Flux_Data},amr::AMR{2,2}) # without face area
    _,_,midpoint,here_data,there_data = unpack(flux_data)
    Δt = amr.global_data.status.Δt
    here_mid = here_vs.midpoint;there_mid = there_vs.midpoint
    # here_weight = here_vs.weight;there_weight = there_vs.weight
    here_vn = here_vs.vn;there_vn = there_vs.vn
    here_df = here_vs.df;there_df = there_vs.df
    here_sdf = here_vs.sdf;there_sdf = there_vs.sdf
    here_ps_mid = here_data.midpoint;there_ps_mid = there_data.midpoint
    # flux = Vector{Float64}(undef,4)
    @inbounds @views begin
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx[i,:],there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
        # flux[1] = sum(here_weight.*here_vn.*df[:,1])+sum(there_weight.*there_vn.*ndf[:,1])
        # flux[2] = sum(here_weight.*here_vn.*df[:,1].*here_mid[:,1])+sum(there_weight.*there_vn.*ndf[:,1].*there_mid[:,1])
        # flux[3] = sum(here_weight.*here_vn.*df[:,1].*here_mid[:,2])+sum(there_weight.*there_vn.*ndf[:,1].*there_mid[:,2])
        # flux[end] = 0.5*(sum(here_weight.*here_vn.*(sum(x->x^2,here_mid,dims = 2).*df[:,1]+df[:,2]))+
        #     sum(there_weight.*there_vn.*(sum(x->x^2,there_mid,dims = 2).*ndf[:,1]+ndf[:,2])))
        here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
        there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    end
    return nothing,[Δt*here_micro,Δt*there_micro]
end