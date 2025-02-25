function face_area(ps_data::AbstractPsData{2}, DIR::Integer)
    return ps_data.ds[FAT[1][DIR]]
end
function face_area(ps_data::AbstractPsData{3}, DIR::Integer)
    return reduce(*, @view(ps_data.ds[FAT[2][DIR]]))
end
function flux!(face::DomainFace,amr::AMR)
    here_vs = make_face_vs(face)
    flux,micro_flux = calc_domain_flux(amr.global_data.config.solver.flux,here_vs,face,amr)
    update_domain_flux!(flux,micro_flux,face,here_vs.heavi)
end
function flux!(face::FullFace,amr::AMR)
    here_vs,there_vs = make_face_vs(face)
    flux,micro_flux = calc_flux(amr.global_data.config.solver.flux,here_vs,there_vs,face,amr)
    update_flux!(flux,micro_flux,face,here_vs.heavi)
end
function flux!(face::HangingFace{DIM,NDF},amr::AMR) where{DIM,NDF}
    rot,direction,midpoint,here_data,there_data = unpack(face)
    here_vs,there_vs = make_face_vs(face)
    for i in eachindex(there_vs)
        flux_data = Flux_Data{HangingFace{DIM,NDF}}(rot,direction,midpoint[i],here_data,there_data[i])
        flux,micro_flux = calc_flux(amr.global_data.config.solver.flux,here_vs,there_vs[i],flux_data,amr)
        update_flux!(flux,micro_flux,flux_data,here_vs.heavi)
    end
end
function flux!(face::BackHangingFace{DIM,NDF},amr::AMR) where{DIM,NDF}
    rot,direction,midpoint,here_data,there_data = unpack(face)
    here_vs,there_vs = make_face_vs(face)
    for i in eachindex(here_vs)
        flux_data = Flux_Data{BackHangingFace{DIM,NDF}}(rot,direction,midpoint[i],here_data[i],there_data)
        flux,micro_flux = calc_flux(amr.global_data.config.solver.flux,here_vs[i],there_vs,flux_data,amr)
        update_flux!(flux,micro_flux,flux_data,here_vs[i].heavi)
    end
end
function update_domain_flux!(flux::AbstractVector,micro_flux::Vector{Matrix{Float64}},face::DomainFace,heavi::Vector{Bool})
    rot,direction,_,_,ps_data = unpack(face)
    area = rot*face_area(ps_data,direction)
    here_micro,there_micro = micro_flux
    ps_data.flux .+= area*flux
    ps_data.vs_data.flux[heavi,:] .+= area*here_micro
    ps_data.vs_data.flux[.!heavi,:] .+= area*there_micro
end
function update_domain_flux!(::Nothing,micro_flux::Vector{Matrix{Float64}},face::DomainFace,heavi::Vector{Bool})
    rot,direction,_,_,ps_data = unpack(face)
    area = rot*face_area(ps_data,direction)
    here_micro,there_micro = micro_flux
    ps_data.vs_data.flux[heavi,:] .+= area*here_micro
    ps_data.vs_data.flux[.!heavi,:] .+= area*there_micro
end
function face_area(::Flux_Data{HangingFace{DIM,NDF}},here_data::AbstractPsData{DIM},direction,rot) where{DIM,NDF}
    face_area(here_data,direction)/2^(DIM-1)*rot
end
function face_area(::T,here_data,direction,rot) where{T<:Union{Flux_Data{BackHangingFace{DIM,NDF}},FullFace{DIM,NDF}} where{DIM,NDF}}
    face_area(here_data,direction)*rot
end
function update_flux!(flux::AbstractVector,micro_flux::Vector{Matrix{Float64}},face::Union{FullFace,Flux_Data},heavi::Vector{Bool})
    rot,direction,_,here_data,there_data = unpack(face)
    area = face_area(face,here_data,direction,rot)
    here_micro,there_micro = micro_flux
    update_macro_flux!(flux*area,here_data,there_data)
    update_micro_flux!(here_micro*area,there_micro*area,here_data,there_data,heavi)
end
function update_flux!(::Nothing,micro_flux::Vector{Matrix{Float64}},face::Union{FullFace,Flux_Data},heavi::Vector{Bool})
    rot,direction,_,here_data,there_data = unpack(face)
    area = face_area(face,here_data,direction,rot)
    here_micro,there_micro = micro_flux
    update_micro_flux!(here_micro*area,there_micro*area,here_data,there_data,heavi)
end
function update_macro_flux!(flux::Vector,here_data::PS_Data,there_data::PS_Data)
    here_data.flux .+= flux
    if there_data.bound_enc>=0
        there_data.flux .-= flux
    end
end
function update_macro_flux!(flux::Vector,here_data::PS_Data,::AbstractGhostPsData)
    here_data.flux .+= flux
end
function update_micro_flux!(here_micro,there_micro,here_data::PS_Data{DIM},::SolidNeighbor{DIM},heavi::Vector{Bool}) where{DIM}
    vs_data = here_data.vs_data;flux = vs_data.flux
    @. flux[heavi,:]+=@views here_micro
    nheavi = [!x for x in heavi]
    @. flux[nheavi,:]+=@views there_micro
end
function update_micro_flux!(here_micro,there_micro,here_data::PS_Data{DIM},there_data::PS_Data,heavi::Vector{Bool}) where{DIM}
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    level = vs_data.level;level_n = nvs_data.level
    flux = vs_data.flux;flux_n = nvs_data.flux
    index = j = index_n = 1;flag = 0.
    if there_data.bound_enc<0
        @inbounds for i = 1:vs_data.vs_num
            if heavi[i]
                @. flux[i, :] += @views here_micro[index, :]
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
                    @. flux[i, :] += @views there_micro[index_n, :]
                    j += 1
                    index_n += 1
                elseif level[i] < level_n[j]
                    while flag != 1.0
                        @. flux[i, :] +=
                            @view(there_micro[index_n, :]) /
                            2^(DIM * (level_n[j] - level[i]))
                        flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                        j += 1
                        index_n += 1
                    end
                    flag = 0.0
                else
                    @. flux[i, :] += @views there_micro[index_n, :]
                    flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                    if flag == 1.0
                        j += 1
                        index_n += 1
                        flag = 0.0
                    end
                end
            end
        end
    else
        @inbounds for i = 1:vs_data.vs_num
            if heavi[i]
                @. flux[i, :] += @views here_micro[index, :]
                if level[i] == level_n[j]
                    @. flux_n[j, :] -= @views here_micro[index, :]
                    j += 1
                elseif level[i] < level_n[j]
                    while flag != 1.0
                        @. flux_n[j, :] -= @views here_micro[index, :]
                        flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                        j += 1
                    end
                    flag = 0.0
                else
                    @. flux_n[j, :] -=
                        @view(here_micro[index, :]) / 2^(DIM * (level[i] - level_n[j]))
                    flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                    if flag == 1.0
                        j += 1
                        flag = 0.0
                    end
                end
                index += 1
            else
                if level[i] == level_n[j]
                    @. flux_n[j, :] -= @views there_micro[index_n, :]
                    @. flux[i, :] += @views there_micro[index_n, :]
                    j += 1
                    index_n += 1
                elseif level[i] < level_n[j]
                    while flag != 1.0
                        @. flux_n[j, :] -= @views there_micro[index_n, :]
                        @. flux[i, :] +=
                            @view(there_micro[index_n, :]) /
                            2^(DIM * (level_n[j] - level[i]))
                        flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                        j += 1
                        index_n += 1
                    end
                    flag = 0.0
                else
                    @. flux[i, :] += @views there_micro[index_n, :]
                    flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                    if flag == 1.0
                        @. flux_n[j, :] -= @views there_micro[index_n, :]
                        j += 1
                        index_n += 1
                        flag = 0.0
                    end
                end
            end
        end
    end
end
function update_micro_flux!(here_micro,there_micro,here_data::PS_Data{DIM},there_data::AbstractGhostPsData,heavi::Vector{Bool}) where{DIM}
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    level = vs_data.level;level_n = nvs_data.level
    flux = vs_data.flux
    index = j = index_n = 1;flag = 0.
    @inbounds for i = 1:vs_data.vs_num
        if heavi[i]
            @. flux[i, :] += @views here_micro[index, :]
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
                @. flux[i, :] += @views there_micro[index_n, :]
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    @. flux[i, :] +=
                        @view(there_micro[index_n, :]) /
                        2^(DIM * (level_n[j] - level[i]))
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                @. flux[i, :] += @views there_micro[index_n, :]
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
end

function make_face_vs(face::DomainFace{DIM,NDF})where{DIM,NDF}
    rot,direction,_,_,ps_data = unpack(face)
    vs_data = ps_data.vs_data;here_mid = vs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    @views here_vs = Face_VS_Data{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],
        vs_data.df[heavi,:],vs_data.sdf[heavi,:,:]
    )
    return here_vs
end
function make_face_vs(face::FullFace{DIM,NDF}) where{DIM,NDF}
    rot,direction,_,here_data,there_data = unpack(face)
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    here_mid = vs_data.midpoint;there_mid = nvs_data.midpoint
    there_data.bound_enc<0 && calc_solid_cell_slope!(nvs_data,vs_data,there_data.midpoint,here_data.midpoint,direction)
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    nheavi = [x>0. for x in rot.*@views there_mid[:,direction]]
    @views here_vs = Face_VS_Data{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],
        vs_data.df[heavi,:],vs_data.sdf[heavi,:,:]
    )
    @views there_vs = Face_VS_Data{DIM,NDF}(
        nheavi,nvs_data.weight[nheavi],there_mid[nheavi,:],there_mid[nheavi,direction],
        nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
    )
    return here_vs,there_vs
end
function make_face_vs(face::HangingFace{DIM,NDF}) where{DIM,NDF}
    rot,direction,_,here_data,there_data = unpack(face)
    vs_data = here_data.vs_data
    here_mid = vs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    @views here_vs = Face_VS_Data{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],vs_data.df[heavi,:],
        vs_data.sdf[heavi,:,:]
    )
    there_vs = Vector{Face_VS_Data{DIM,NDF}}(undef,length(there_data))
    for i in eachindex(there_data)
        nvs_data = there_data[i].vs_data
        there_mid = nvs_data.midpoint
        nheavi = [x>0. for x in rot.*@views there_mid[:,direction]]
        @views there_vs[i] = Face_VS_Data{DIM,NDF}(
            nheavi,nvs_data.weight[nheavi],there_mid[nheavi,:],there_mid[nheavi,direction],
            nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
        )
    end
    return here_vs,there_vs
end
function make_face_vs(face::BackHangingFace{DIM,NDF}) where{DIM,NDF}
    rot,direction,_,here_data,there_data = unpack(face)
    nvs_data = there_data.vs_data
    nheavi = [x>0. for x in rot.*@views nvs_data.midpoint[:,direction]]
    @views there_vs = Face_VS_Data{DIM,NDF}(
            nheavi,nvs_data.weight[nheavi],nvs_data.midpoint[nheavi,:],nvs_data.midpoint[nheavi,direction],
            nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
        )
    here_vs = Vector{Face_VS_Data{DIM,NDF}}(undef,length(here_data))
    for i in eachindex(here_data)
        vs_data = here_data[i].vs_data
        here_mid = vs_data.midpoint
        heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
        @views here_vs[i] = Face_VS_Data{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],vs_data.df[heavi,:],
        vs_data.sdf[heavi,:,:]
    )
    end
    return here_vs,there_vs
end
function flux!(amr::AMR)
    faces = amr.field.faces
    @simd for face in faces
        flux!(face,amr)
    end
end