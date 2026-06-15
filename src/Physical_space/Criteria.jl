const PS_COARSEN_SENSOR_RATIO = 0.3
const PS_LOHNER_ABS_FLOOR = 1e-4
const PS_PRIMITIVE_REL_JUMP_FLOOR = 1e-3
const PS_VORTICITY_JUMP_FLOOR = 2e-2
const PS_BOUNDARY_SENSOR_FRACTIONS = (-0.5, -0.25, 0.0, 0.25, 0.5)

@inline function lohner_rowmax(lohner::AbstractMatrix, row::Integer)
    value = zero(eltype(lohner))
    @inbounds for dir in axes(lohner, 2)
        value = max(value, lohner[row, dir])
    end
    return value
end

@inline function ps_sensor(ps_data::PsData)
    max(
        lohner_rowmax(ps_data.lohner, 1),
        # max(
            # lohner_rowmax(ps_data.lohner, 2),
            lohner_rowmax(ps_data.lohner, size(ps_data.lohner, 1))
        # ),
    )
end

@inline function lohner_value(left::Real, center::Real, right::Real, dsL::Real, dsR::Real, eps::Real)
    scale = dsR * abs(left) + (dsL + dsR) * abs(center) + dsL * abs(right)
    scale <  PS_LOHNER_ABS_FLOOR*min(dsL,dsR) && return 0.0
    denom = dsR * abs(left - center) + dsL * abs(right - center) + eps * scale
    denom <= 0.0 && return 0.0
    abs(dsR * left - (dsL + dsR) * center + dsL * right) / denom
end

@inline function primitive_amplitude_ok(left::Real, center::Real, right::Real)
    jump = max(abs(left - center), abs(right - center))
    scale = max(abs(center), PS_LOHNER_ABS_FLOOR)
    jump >= PS_PRIMITIVE_REL_JUMP_FLOOR * scale
end

@inline function velocity_slope(sw::AbstractMatrix, prim::AbstractVector, component::Integer, dir::Integer)
    row = component + 1
    @inbounds (sw[row, dir] - prim[row] * sw[1, dir]) / prim[1]
end

@inline function vorticity(sw::AbstractMatrix, prim::AbstractVector, ::Val{2})
    velocity_slope(sw, prim, 1, 2) - velocity_slope(sw, prim, 2, 1)
end

@inline function vorticity(sw::AbstractMatrix, prim::AbstractVector, ::Val{3})
    c1 = velocity_slope(sw, prim, 2, 3) - velocity_slope(sw, prim, 3, 2)
    c2 = velocity_slope(sw, prim, 3, 1) - velocity_slope(sw, prim, 1, 3)
    c3 = velocity_slope(sw, prim, 1, 2) - velocity_slope(sw, prim, 2, 1)
    sqrt(c1 * c1 + c2 * c2 + c3 * c3)
end

@inline function velocity_scale(prim::AbstractVector)
    speed2 = 0.0
    @inbounds for i in 2:(length(prim) - 1)
        speed2 += prim[i] * prim[i]
    end
    lambda = max(abs(prim[end]), eps(Float64))
    max(sqrt(speed2), inv(sqrt(lambda)))
end

@inline function vorticity_amplitude_ok(left::Real, center::Real, right::Real, prim::AbstractVector, h::Real)
    omega = max(abs(left), max(abs(center), abs(right)))
    omega * h >= PS_VORTICITY_JUMP_FLOOR * velocity_scale(prim)
end

@inline _ps_domain_face_id(dir::Integer, boundary_on_left::Bool) =
    boundary_on_left ? 2 * dir - 1 : 2 * dir

function _ps_domain_for_face(kinfo::KInfo, faceid::Integer)
    domains = kinfo.config.domain
    faceid <= length(domains) && domains[faceid].id == faceid && return domains[faceid]
    i = findfirst(domain -> domain.id == faceid, domains)
    return i === nothing ? nothing : domains[i]
end

function _ps_domain_bc_prim(domain::Domain, midpoint::AbstractVector, kinfo::KInfo)
    isdefined(domain, :bc) || return nothing
    bc = domain.bc
    if bc isa AbstractVector
        return Float64.(collect(bc))
    elseif bc isa Function
        if applicable(bc, midpoint, kinfo)
            return Float64.(collect(bc(midpoint, kinfo)))
        elseif applicable(bc, midpoint)
            return Float64.(collect(bc(midpoint)))
        end
        return Float64.(collect(bc(; midpoint)))
    end
    return nothing
end

function _ps_boundary_sample_point!(
    face_midpoint::AbstractVector,
    ps_data::PsData{DIM,NDF},
    kinfo::KInfo,
    faceid::Integer,
    dir::Integer,
    sample::Integer,
) where {DIM,NDF}
    face_midpoint .= ps_data.midpoint
    face_midpoint[dir] = kinfo.config.geometry[faceid]

    q = sample - 1
    nfrac = length(PS_BOUNDARY_SENSOR_FRACTIONS)
    @inbounds for d in 1:DIM
        d == dir && continue
        i = mod(q, nfrac) + 1
        q = div(q, nfrac)
        x = ps_data.midpoint[d] + PS_BOUNDARY_SENSOR_FRACTIONS[i] * ps_data.ds[d]
        face_midpoint[d] = clamp(x, kinfo.config.geometry[2d - 1], kinfo.config.geometry[2d])
    end
    return face_midpoint
end

@inline function _ps_valid_neighbor(data)
    return data !== nothing && !isa(data, InsideSolidData) && !isa(data, GhostInsideSolidData)
end

function _ps_average_neighbor_prim!(
    prim::Vector{Float64},
    neighbor_data::AbstractVector,
    ps_data::PsData{DIM,NDF},
    dir::Integer,
    ds_neighbor::Real,
    kinfo::KInfo,
) where {DIM,NDF}
    fill!(prim, 0.0)
    n = 0
    first_neighbor = nothing
    @inbounds for data in neighbor_data
        _ps_valid_neighbor(data) || continue
        first_neighbor === nothing && (first_neighbor = data)
        @. prim += data.w
        n += 1
    end
    if n == 0
        prim .= ps_data.prim
        return prim
    end

    if ds_neighbor > ps_data.ds[dir] && first_neighbor !== nothing
        dx = (ps_data.midpoint - first_neighbor.midpoint)[FAT[DIM - 1][dir]]
        @inbounds for j in eachindex(prim)
            @views prim[j] += dot(dx, first_neighbor.sw[j, FAT[DIM - 1][dir]])
        end
    end

    prim ./= n
    prim .= get_prim(prim, kinfo)
    return prim
end

function update_Lohner_boundary_ps!(
    ps_data::PsData{DIM,NDF},
    interior_data::AbstractVector,
    ds_boundary::Real,
    ds_interior::Real,
    dir::Integer,
    boundary_on_left::Bool,
    ws_boundary::Vector{Float64},
    ws_interior::Vector{Float64},
    kinfo::KInfo,
) where {DIM,NDF}
    faceid = _ps_domain_face_id(dir, boundary_on_left)
    domain = _ps_domain_for_face(kinfo, faceid)
    if domain === nothing || !isdefined(domain, :bc)
        ps_data.lohner[:, dir] .= 0.0
        return nothing
    end

    _ps_average_neighbor_prim!(ws_interior, interior_data, ps_data, dir, ds_interior, kinfo)
    dsL = boundary_on_left ? ds_boundary : ds_interior
    dsR = boundary_on_left ? ds_interior : ds_boundary
    eps_l = 0.2 * ps_data.ds[dir]

    ps_data.lohner[:, dir] .= 0.0
    face_midpoint = similar(ps_data.midpoint)
    nsamples = length(PS_BOUNDARY_SENSOR_FRACTIONS)^(DIM - 1)
    @inbounds for sample in 1:nsamples
        _ps_boundary_sample_point!(face_midpoint, ps_data, kinfo, faceid, dir, sample)
        boundary_prim = _ps_domain_bc_prim(domain, face_midpoint, kinfo)
        boundary_prim === nothing && continue
        length(boundary_prim) == length(ps_data.prim) ||
            error("Domain boundary $faceid returned a primitive vector of length $(length(boundary_prim)); expected $(length(ps_data.prim)).")
        ws_boundary .= boundary_prim

        left = boundary_on_left ? ws_boundary : ws_interior
        right = boundary_on_left ? ws_interior : ws_boundary
        for j in eachindex(ps_data.prim)
            value = primitive_amplitude_ok(left[j], ps_data.prim[j], right[j]) ?
                lohner_value(left[j], ps_data.prim[j], right[j], dsL, dsR, eps_l) : 0.0
            ps_data.lohner[j, dir] = max(ps_data.lohner[j, dir], value)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function ps_refine_flag(
    ps_data::PsData{DIM},
    level::Int8,
    ka::KA{DIM}
) where{DIM}
    kinfo = ka.kinfo
    if ps_data.bound_enc!=0||domain_flag(kinfo,ps_data.midpoint,ps_data.ds)
        return Cint(1)
    end
    level>kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL-1&&return Cint(0)
    kinfo.config.user_defined.static_ps_refine_flag(ps_data.midpoint,ps_data.ds,kinfo,level) && return Cint(1)
    dflag = kinfo.config.user_defined.dynamic_ps_refine_flag==null_udf ? true : kinfo.config.user_defined.dynamic_ps_refine_flag(ps_data,level,ka)
    !dflag&&return Cint(0)
    return Cint(ps_sensor(ps_data)>kinfo.config.solver.ADAPT_COEFFI_PS)
end

"""
$(TYPEDSIGNATURES)
"""
function ps_coarsen_flag(ps_datas::Vector{PsData}, levels::Vector{Int}, ka::KA{DIM,NDF}) where{DIM,NDF}
    kinfo = ka.kinfo
    levels[1]>kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL&&return Cint(0)
    threshold = PS_COARSEN_SENSOR_RATIO * kinfo.config.solver.ADAPT_COEFFI_PS
    for i = 1:2^DIM
        ps_data = ps_datas[i]
        (ps_data.bound_enc!=0||domain_flag(kinfo,ps_data.midpoint,ps_data.ds)) && return Cint(0)
        kinfo.config.user_defined.static_ps_refine_flag(ps_data.midpoint,ps_data.ds,kinfo,levels[i]-1) && return Cint(0)
        ps_sensor(ps_data)>threshold && return Cint(0)
    end
    return Cint(1)
end

function update_Lohner_inner_ps!(
    ps_data::PsData{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
    kinfo::KInfo
) where{DIM,NDF}
    fill!(ws_swL, 0.0)
    fill!(ws_swR, 0.0)
    swL = zeros(DIM+2,DIM)
    swR = zeros(DIM+2,DIM)
    nL = length(Ldata)
    nR = length(Rdata)

    @inbounds for j in 1:nL
        @. ws_swL += Ldata[j].w
        @. swL += Ldata[j].sw
    end
    @inbounds for j in 1:nR
        @. ws_swR += Rdata[j].w
        @. swR += Rdata[j].sw
    end

    if dsL>ps_data.ds[dir]
        dx = (ps_data.midpoint-Ldata[1].midpoint)[FAT[DIM-1][dir]]
        @inbounds for j in eachindex(ws_swL)
            @views ws_swL[j] += dot(dx, Ldata[1].sw[j,FAT[DIM-1][dir]])
        end
    end
    if dsR>ps_data.ds[dir]
        dx = (ps_data.midpoint-Rdata[1].midpoint)[FAT[DIM-1][dir]]
        @inbounds for j in eachindex(ws_swR)
            @views ws_swR[j] += dot(dx, Rdata[1].sw[j,FAT[DIM-1][dir]])
        end
    end

    ws_swL ./= nL
    ws_swR ./= nR
    swL ./= nL
    swR ./= nR

    ws_swL .= get_prim(ws_swL,kinfo)
    ws_swR .= get_prim(ws_swR,kinfo)

    omegaL = vorticity(swL, ws_swL, Val(DIM))
    omegaR = vorticity(swR, ws_swR, Val(DIM))
    omega = vorticity(ps_data.sw, ps_data.prim, Val(DIM))
    use_vorticity = vorticity_amplitude_ok(omegaL, omega, omegaR, ps_data.prim, max(dsL, dsR))

    eps_l = 0.2*ps_data.ds[dir]
    @inbounds for j in eachindex(ws_swL)
        if j == 2
            ps_data.lohner[j, dir] = use_vorticity ?
                lohner_value(omegaL, omega, omegaR, dsL, dsR, eps_l) : 0.0
        else
            ps_data.lohner[j, dir] =
                primitive_amplitude_ok(ws_swL[j], ps_data.prim[j], ws_swR[j]) ?
                lohner_value(ws_swL[j], ps_data.prim[j], ws_swR[j], dsL, dsR, eps_l) : 0.0
        end
    end

    return nothing
end
