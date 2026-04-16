include("CAIDVM.jl")
include("DVM.jl")
include("UGKS.jl")
include("Slope.jl")
export update_slope!, update_slope_inner_vs!, update_slope_bound_vs!, update_slope_inner_ps!, update_slope_bound_ps!
export vanleer, minmod, diff_vs!
export flux!, update_flux!, update_domain_flux!, update_micro_flux!, update_macro_flux!, make_face_vs
export calc_flux, calc_domain_flux, positivity_preserving_reconstruct

function face_area(ps_data::AbstractPsData{2}, DIR::Integer)
    return ps_data.ds[FAT[1][DIR]]
end
function face_area(ps_data::AbstractPsData{3}, DIR::Integer)
    return reduce(*, @view(ps_data.ds[FAT[2][DIR]]))
end
"""
$(TYPEDSIGNATURES)
"""
function flux!(face::DomainFace,ka::KA)
    here_vs = make_face_vs(face)
    flux,micro_flux = calc_domain_flux(ka.kinfo.config.solver.flux,here_vs,face,ka)
    update_domain_flux!(flux,micro_flux,face,here_vs.heavi)
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function flux!(face::FullFace,ka::KA)
    here_vs,there_vs = make_face_vs(face)
    flux,micro_flux = calc_flux(ka.kinfo.config.solver.flux,here_vs,there_vs,face,ka)
    update_flux!(flux,micro_flux,face,here_vs.heavi)
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function flux!(face::HangingFace{DIM,NDF},ka::KA) where{DIM,NDF}
    rot,direction,midpoint,here_data,there_data = unpack(face)
    here_vs,there_vs = make_face_vs(face)
    for i in eachindex(there_vs)
        flux_data = FluxData{HangingFace{DIM,NDF}}(rot,direction,midpoint[i],here_data,there_data[i])
        flux,micro_flux = calc_flux(ka.kinfo.config.solver.flux,here_vs,there_vs[i],flux_data,ka)
        update_flux!(flux,micro_flux,flux_data,here_vs.heavi)
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function flux!(face::BackHangingFace{DIM,NDF},ka::KA) where{DIM,NDF}
    rot,direction,midpoint,here_data,there_data = unpack(face)
    here_vs,there_vs = make_face_vs(face)
    for i in eachindex(here_vs)
        flux_data = FluxData{BackHangingFace{DIM,NDF}}(rot,direction,midpoint[i],here_data[i],there_data)
        flux,micro_flux = calc_flux(ka.kinfo.config.solver.flux,here_vs[i],there_vs,flux_data,ka)
        update_flux!(flux,micro_flux,flux_data,here_vs[i].heavi)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_domain_flux!(flux::AbstractVector,micro_flux::Vector{Matrix{Float64}},face::DomainFace,heavi::Vector{Bool})
    rot,direction,_,_,ps_data = unpack(face)
    area = rot*face_area(ps_data,direction)
    here_micro,there_micro = micro_flux
    ps_data.flux .+= area*flux
    ps_data.vs_data.flux[heavi,:] .+= area*here_micro
    ps_data.vs_data.flux[.!heavi,:] .+= area*there_micro
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_domain_flux!(::Nothing,micro_flux::Vector{Matrix{Float64}},face::DomainFace,heavi::Vector{Bool})
    rot,direction,_,_,ps_data = unpack(face)
    area = rot*face_area(ps_data,direction)
    here_micro,there_micro = micro_flux
    ps_data.vs_data.flux[heavi,:] .+= area*here_micro
    ps_data.vs_data.flux[.!heavi,:] .+= area*there_micro
    return nothing
end
function face_area(::FluxData{HangingFace{DIM,NDF}},here_data::AbstractPsData{DIM},direction,rot) where{DIM,NDF}
    face_area(here_data,direction)/2^(DIM-1)*rot
end
function face_area(::T,here_data,direction,rot) where{T<:Union{FluxData{BackHangingFace{DIM,NDF}},FullFace{DIM,NDF}} where{DIM,NDF}}
    face_area(here_data,direction)*rot
end

"""
$(TYPEDSIGNATURES)
"""
function update_flux!(flux::AbstractVector,micro_flux::Vector{Matrix{Float64}},face::Union{FullFace,FluxData},heavi::Vector{Bool})
    rot,direction,_,here_data,there_data = unpack(face)
    area = face_area(face,here_data,direction,rot)
    here_micro,there_micro = micro_flux
    update_macro_flux!(flux*area,here_data,there_data)
    update_micro_flux!(here_micro*area,there_micro*area,here_data,there_data,heavi)
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_flux!(::Nothing,micro_flux::Vector{Matrix{Float64}},face::Union{FullFace,FluxData},heavi::Vector{Bool})
    rot,direction,_,here_data,there_data = unpack(face)
    area = face_area(face,here_data,direction,rot)
    here_micro,there_micro = micro_flux
    update_micro_flux!(here_micro*area,there_micro*area,here_data,there_data,heavi)
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_macro_flux!(flux::Vector,here_data::PsData,there_data::PsData)
    here_data.flux .+= flux
    if there_data.bound_enc>=0
        there_data.flux .-= flux
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_macro_flux!(flux::Vector,here_data::PsData,::AbstractGhostPsData)
    here_data.flux .+= flux
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_macro_flux!(flux::Vector,here_data::PsData,::SolidNeighbor)
    here_data.flux .+= flux
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_micro_flux!(here_micro,there_micro,here_data::PsData{DIM},::SolidNeighbor{DIM},heavi::Vector{Bool}) where{DIM}
    vs_data = here_data.vs_data;flux = vs_data.flux
    @. flux[heavi,:]+= here_micro
    nheavi = [!x for x in heavi]
    @. flux[nheavi,:]+= there_micro
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_micro_flux!(here_micro,there_micro,here_data::PsData{DIM,NDF},there_data::PsData,heavi::Vector{Bool}) where{DIM,NDF}
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    level = vs_data.level;level_n = nvs_data.level
    flux = vs_data.flux;flux_n = nvs_data.flux
    index = j = index_n = 1;flag = 0.
    @inbounds for i = 1:vs_data.vs_num
        if heavi[i]
            @simd for ii in 1:NDF
                flux[i, ii] +=here_micro[index, ii]
            end
            if level[i] == level_n[j]
                @simd for ii in 1:NDF
                    flux_n[j, ii] -= here_micro[index, ii]
                end
                j += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    @simd for ii in 1:NDF
                        flux_n[j, ii] -= here_micro[index, ii]
                    end
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                end
                flag = 0.0
            else
                @simd for ii in 1:NDF
                    flux_n[j, ii] -=
                        (here_micro[index, ii]) / 2^(DIM * (level[i] - level_n[j]))
                end
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
            index += 1
        else
            if level[i] == level_n[j]
                @simd for ii in 1:NDF
                    flux_n[j, ii] -= there_micro[index_n, ii]
                    flux[i, ii] += there_micro[index_n, ii]
                end
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    @simd for ii in 1:NDF
                        flux_n[j, ii] -= there_micro[index_n, ii]
                        flux[i, ii] +=
                            (there_micro[index_n, ii]) /
                            2^(DIM * (level_n[j] - level[i]))
                    end
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                @simd for ii in 1:NDF
                    flux[i, ii] +=  there_micro[index_n, ii]
                end
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    @simd for ii in 1:NDF
                        flux_n[j, ii] -= there_micro[index_n, ii]
                    end
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_micro_flux!(here_micro,there_micro,here_data::PsData{DIM,NDF},there_data::AbstractGhostPsData,heavi::Vector{Bool}) where{DIM,NDF}
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    level = vs_data.level;level_n = nvs_data.level
    flux = vs_data.flux
    index = j = index_n = 1;flag = 0.
    @inbounds for i = 1:vs_data.vs_num
        if heavi[i]
            @simd for ii in 1:NDF
                flux[i, ii] +=here_micro[index, ii]
            end
            if level[i] == level_n[j]
                j += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                end
                flag = 0.0
            else
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
            index += 1
        else
            if level[i] == level_n[j]
                @simd for ii in 1:NDF
                    flux[i, ii] += there_micro[index_n, ii]
                end
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    for ii in 1:NDF
                        flux[i, ii] +=
                            (there_micro[index_n, ii]) /
                            2^(DIM * (level_n[j] - level[i]))
                    end
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                for ii in 1:NDF
                    flux[i, ii] += there_micro[index_n, ii]
                end
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function make_face_vs(face::DomainFace{DIM,NDF})where{DIM,NDF}
    rot,direction,_,_,ps_data = unpack(face)
    vs_data = ps_data.vs_data;here_mid = vs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    @views here_vs = FaceVsData{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],
        vs_data.df[heavi,:],vs_data.sdf[heavi,:,:]
    )
    return here_vs
end
"""
$(TYPEDSIGNATURES)
"""
function make_face_vs(face::FullFace{DIM,NDF}) where{DIM,NDF}
    rot,direction,_,here_data,there_data = unpack(face)
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    here_mid = vs_data.midpoint;there_mid = nvs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    nheavi = [x>0. for x in rot.*@views there_mid[:,direction]]
    @views here_vs = FaceVsData{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],
        vs_data.df[heavi,:],vs_data.sdf[heavi,:,:]
    )
    @views there_vs = FaceVsData{DIM,NDF}(
        nheavi,nvs_data.weight[nheavi],there_mid[nheavi,:],there_mid[nheavi,direction],
        nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
    )
    return here_vs,there_vs
end
"""
$(TYPEDSIGNATURES)
"""
function make_face_vs(face::HangingFace{DIM,NDF}) where{DIM,NDF}
    rot,direction,_,here_data,there_data = unpack(face)
    vs_data = here_data.vs_data
    here_mid = vs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    @views here_vs = FaceVsData{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],vs_data.df[heavi,:],
        vs_data.sdf[heavi,:,:]
    )
    there_vs = Vector{FaceVsData{DIM,NDF}}(undef,length(there_data))
    for i in eachindex(there_data)
        nvs_data = there_data[i].vs_data
        there_mid = nvs_data.midpoint
        nheavi = [x>0. for x in rot.*@views there_mid[:,direction]]
        @views there_vs[i] = FaceVsData{DIM,NDF}(
            nheavi,nvs_data.weight[nheavi],there_mid[nheavi,:],there_mid[nheavi,direction],
            nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
        )
    end
    return here_vs,there_vs
end
"""
$(TYPEDSIGNATURES)
"""
function make_face_vs(face::BackHangingFace{DIM,NDF}) where{DIM,NDF}
    rot,direction,_,here_data,there_data = unpack(face)
    nvs_data = there_data.vs_data
    nheavi = [x>0. for x in rot.*@views nvs_data.midpoint[:,direction]]
    @views there_vs = FaceVsData{DIM,NDF}(
            nheavi,nvs_data.weight[nheavi],nvs_data.midpoint[nheavi,:],nvs_data.midpoint[nheavi,direction],
            nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
        )
    here_vs = Vector{FaceVsData{DIM,NDF}}(undef,length(here_data))
    for i in eachindex(here_data)
        vs_data = here_data[i].vs_data
        here_mid = vs_data.midpoint
        heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
        @views here_vs[i] = FaceVsData{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],vs_data.df[heavi,:],
        vs_data.sdf[heavi,:,:]
    )
    end
    return here_vs,there_vs
end

"""
$(TYPEDSIGNATURES)
Outer function for flux computation and update. The iteration is carried out through faces, which avoids redundant computation.
"""
function flux!(ka::KA)
    faces = ka.kdata.field.faces
    for face in faces
        flux!(face,ka)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Asynchronous flux computation with overlapped `solid_exchange!`.

1. `update_solid_cell!` — interpolate solid-cell df from fluid neighbors.
2. `solid_exchange_begin!` — post non-blocking MPI sends/recvs.
3. Compute flux for all non-IB faces (overlapped with communication).
4. `solid_exchange_finish!` — `MPI_Waitall` + scatter.
5. `update_solid_neighbor!` — reconstruct `SolidNeighbor` df/w/sw from the
   now-synchronised ghost data.
6. Compute flux for IB faces.
"""
function flux!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where {DIM,NDF}
    ibs = ka.kdata.field.immersed_boundaries

    # 1. update solid cells (local interpolation)
    update_solid_cell!(ka)

    # 2. begin async solid exchange
    recv_bufs = solid_exchange_begin!(p4est, ka)

    # 3. compute non-IB face fluxes (overlapped with MPI communication)
    ib_face_set = Set{UInt}()
    for ib in ibs
        for f in ib.faces; push!(ib_face_set, objectid(f)); end
    end
    for face in ka.kdata.field.faces
        objectid(face) in ib_face_set && continue
        flux!(face, ka)
    end

    # 4. finish async solid exchange (Waitall + scatter)
    solid_exchange_finish!(ka, recv_bufs)

    # 5. update solid neighbors
    update_solid_neighbor!(ka)

    # 6. compute IB face fluxes
    for ib in ibs
        for face in ib.faces
            flux!(face, ka)
        end
    end

    return nothing
end