# Helper: reconstruct df from here_df + slope correction without materialising dx.
# df[i,j] += Σ_k (midpoint[k] - ps_mid[k] - mid[i,k]*Δt) * sdf[i,j,k]
@inline function _reconstruct!(df, mid, sdf, ps_mid, midpoint, Δt)
    for k in axes(sdf, 3)
        δxk = midpoint[k] - ps_mid[k]
        @views @. df += (δxk - mid[:, k] * Δt) * sdf[:, :, k]
    end
end

function calc_domain_flux(::DVM,here_vs::FaceVsData,face::DomainFace{DIM,NDF,Maxwellian},ka::KA) where{DIM,NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,here_weight,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = ka.kinfo.status.Δt
    vs_data = ps_data.vs_data
    nheavi = .!heavi
    there_mid    = @views vs_data.midpoint[nheavi,:]
    there_weight = @views vs_data.weight[nheavi]
    there_vn     = @views there_mid[:,direction]
    df = copy(here_df)
    _reconstruct!(df, here_mid, here_sdf, ps_data.midpoint, midpoint, Δt)
    bc = get_bc(face.domain.bc)
    bc[1] = calc_ρw(there_mid,df,bc,here_vn,there_vn,here_weight,there_weight,ka)
    there_df = Matrix{Float64}(undef,size(there_mid,1),NDF)
    for i in axes(there_df,1)
        @views there_df[i,:] .= discrete_maxwell(there_mid[i,:],bc,ka.kinfo)
    end
    here_micro  = df .* here_vn
    there_micro = there_df .* there_vn
    return nothing,[here_micro,there_micro]
end
function calc_domain_flux(::DVM,here_vs::FaceVsData,face::DomainFace{DIM,NDF,SuperSonicInflow},ka::KA) where{DIM,NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,_,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = ka.kinfo.status.Δt
    vs_data = ps_data.vs_data
    nheavi = .!heavi
    there_mid = @views vs_data.midpoint[nheavi,:]
    there_vn  = @views there_mid[:,direction]
    df = copy(here_df)
    _reconstruct!(df, here_mid, here_sdf, ps_data.midpoint, midpoint, Δt)
    bc = get_bc(face.domain.bc)
    there_df = Matrix{Float64}(undef,size(there_mid,1),NDF)
    for i in axes(there_df,1)
        @views there_df[i,:] .= discrete_maxwell(there_mid[i,:],bc,ka.kinfo)
    end
    here_micro  = df .* here_vn
    there_micro = there_df .* there_vn
    return nothing,[here_micro,there_micro]
end
function calc_domain_flux(::DVM,here_vs::FaceVsData,face::DomainFace{DIM,NDF,UniformOutflow},ka::KA) where{DIM,NDF}
    _,direction,_,_,ps_data = unpack(face)
    heavi,_,_,here_vn,_,_ = unpack(here_vs)
    vs_data = ps_data.vs_data
    nheavi = .!heavi
    there_mid = @views vs_data.midpoint[nheavi,:]
    ndf       = @views vs_data.df[nheavi,:]
    there_vn  = @views there_mid[:,direction]
    df = @views vs_data.df[heavi,:]
    here_micro  = df .* here_vn
    there_micro = ndf .* there_vn
    return nothing,[here_micro,there_micro]
end
function calc_domain_flux(::DVM,here_vs::FaceVsData,face::DomainFace{2,NDF,InterpolatedOutflow},ka::KA) where{NDF}
    _,direction,midpoint,_,ps_data = unpack(face)
    heavi,_,here_mid,here_vn,here_df,here_sdf = unpack(here_vs)
    Δt = ka.kinfo.status.Δt
    here_ps_mid  = ps_data.midpoint
    there_ps_mid = 2.0*midpoint - here_ps_mid
    vs_data = ps_data.vs_data
    nheavi = .!heavi
    there_mid = @views vs_data.midpoint[nheavi,:]
    there_sdf = @views vs_data.sdf[nheavi,:,:]
    there_vn  = @views there_mid[:,direction]
    # Extrapolate to the ghost PS midpoint, then apply slope correction.
    there_df0 = vs_data.df[nheavi,:] .+
                (there_ps_mid[direction] - here_ps_mid[direction]) .*
                vs_data.sdf[nheavi,:,direction]
    df  = copy(here_df)
    ndf = there_df0           # already allocated above
    _reconstruct!(df,  here_mid,  here_sdf,  here_ps_mid,  midpoint, Δt)
    _reconstruct!(ndf, there_mid, there_sdf, there_ps_mid, midpoint, Δt)
    here_micro  = df  .* here_vn
    there_micro = ndf .* there_vn
    return nothing,[here_micro,there_micro]
end
function calc_flux(::DVM,here_vs,there_vs,flux_data::Union{FullFace,FluxData},ka::KA{DIM,NDF}) where{DIM,NDF}
    _,_,midpoint,here_data,there_data = unpack(flux_data)
    Δt = ka.kinfo.status.Δt
    here_mid    = here_vs.midpoint;  there_mid    = there_vs.midpoint
    here_vn     = here_vs.vn;        there_vn     = there_vs.vn
    here_df     = here_vs.df;        there_df     = there_vs.df
    here_sdf    = here_vs.sdf;       there_sdf    = there_vs.sdf
    here_ps_mid = here_data.midpoint; there_ps_mid = there_data.midpoint
    df  = copy(here_df)
    ndf = copy(there_df)
    _reconstruct!(df,  here_mid,  here_sdf,  here_ps_mid,  midpoint, Δt)
    _reconstruct!(ndf, there_mid, there_sdf, there_ps_mid, midpoint, Δt)
    here_micro  = df  .* here_vn
    there_micro = ndf .* there_vn
    return nothing,[here_micro,there_micro]
end
