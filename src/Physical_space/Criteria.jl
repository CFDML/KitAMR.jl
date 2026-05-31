const PS_COARSEN_SENSOR_RATIO = 0.3
const PS_LOHNER_ABS_FLOOR = 1e-4
const PS_PRIMITIVE_REL_JUMP_FLOOR = 1e-3
const PS_VORTICITY_JUMP_FLOOR = 2e-2

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
        max(
            lohner_rowmax(ps_data.lohner, 2),
            lohner_rowmax(ps_data.lohner, size(ps_data.lohner, 1)),
        ),
    )
end

@inline function lohner_value(left::Real, center::Real, right::Real, dsL::Real, dsR::Real, eps::Real)
    scale = dsR * abs(left) + (dsL + dsR) * abs(center) + dsL * abs(right)
    scale < PS_LOHNER_ABS_FLOOR && return 0.0
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
