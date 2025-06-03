function update_micro_flux!(here_micro,there_micro,here_data::PS_Data{DIM,NDF},there_data::AbstractGhostPsData,heavi::Vector{Bool}) where{DIM,NDF}
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
end

function solid_cell_exchange!(ps4est,amr)
    solic_cell_exchange!(amr.global_data.config.flux,ps4est,amr)
end
function solid_cell_exchange!(::Type{MicroFlux},ps4est,amr)
    transfer_solid_cell_micro()
end
# function velocity_average_gradient(midpoints::AbstractMatrix,dfs::AbstractMatrix,id::Int,ids::AbstractVector) # Not the most efficient strategy. Maybe abandon the average, use one point for each direction?
#     DIM = size(midpoints,2);NDF = size(dfs,2)
#     vg = Vector{Float64}(undef,NDF,DIM)
#     for i in 1:DIM
#         for j in eachindex(ids)
#             for k in 1:NDF
#             end
#         end
#     end
# end
function velocity_gradient(midpoints::AbstractMatrix,dfs::AbstractMatrix,id::Int,ids::AbstractVector) # Use one point for each direction to achieve higher efficiency. The accuracy is satisfactory because the template is ideally close enough to midpoint_n.
    DIM = size(midpoints,2);NDF = size(dfs,2)
    vg = Vector{Float64}(undef,NDF,DIM)
    for i in 1:DIM
        for j in 1:NDF
            vg[j,i] = (dfs[ids[i],j]-dfs[id,j])/(midpoints[ids[i],i]-midpoints[id,i])
        end
    end
    return vg
end
function push_vg!(heavi_id::Vector{Int},template::Matrix{Int},ps_data::PS_Data,id::Int)
    push!(ps_data.neighbor.vg_template[id],VelocityGradient(template,heavi_id))
    return nothing
end
function push_vg!(::Nothing,::Nothing,::AbstractPsData,::Int)
    return nothing
end
function vg_heavi_id_encode(vg,nvg,vs_data::VS_Data{DIM},nvs_data::VS_Data,direction::Int,rot::Float64) where{DIM}# Should be merged into gradient_encode?
    heavi = [x<=0. for x in rot.*@views vs_data.midpoint[:,direction]]
    level = vs_data.level;level_n = nvs_data.level
    table = zeros(Int,vs_data.vs_num);table_n = zeros(Int,nvs_data.vs_num)
    table_index = table_index_n = 1;index = 1;flag = 0.
    for i in 1:vs_data.vs_num
        if level[i]==level_n[index]
            index+=1
            # root_flag += 1/2^(DIM*level[i])
            if heavi[i]
                table[i] = table_index
                table_index += 1
            else
                table_n[index] = table_index_n
                table_index_n += 1
            end
        elseif level[i]<level_n[index]
            # push!(vg_index,i)
            while flag != 1.0
                # append!(temp_id,search_velocity_gradient_temp(vs_data,i,@views nvs_data.midpoint[index,:],root,nums,global_data))
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                if !heavi[i]
                    table_n[index] = table_index_n
                    table_index_n += 1
                end
                index+=1
            end
            if heavi[i]
                table[i] = table_index
                table_index += 1
            end
            flag = 0.0
        else
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if heavi[i]
                table[i] = table_index
                table_index += 1
            end
            if flag == 1.0
                if !heavi[i] 
                    table_n[index] = table_index_n
                    table_index_n += 1
                end
                index += 1
                flag = 0.0
            end
        end
    end
    index = 1;flag = 0.
    vg_count = vg_count_n = 1
    heavi_index = Int[];heavi_index_n = Int[]
    for i in 1:vs_data.vs_num
        if level[i]==level_n[index]
            index+=1
            # root_flag += 1/2^(DIM*level[i])
        elseif level[i]<level_n[index]
            # push!(vg_index,i)
            while flag != 1.0
                # append!(temp_id,search_velocity_gradient_temp(vs_data,i,@views nvs_data.midpoint[index,:],root,nums,global_data))
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index+=1
            end
            if heavi[i]
                @views append!(heavi_index,table[vg[:,vg_count]])
            end
            vg_count += 1
            flag = 0.0
        else
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                if !heavi[i] 
                    push!(heavi_index_n,table[nvg[:,vg_count_n]])
                end
                vg_count_n += 1
                index += 1
                flag = 0.0
            end
        end
    end
    heavi_index = reshape(heavi_index,DIM,:);heavi_index_n = reshape(heavi_index_n,DIM,:)
    return heavi_index,heavi_index_n
end
function vg_heavi_id_encode(::Nothing,nvg::Matrix,vs_data::VS_Data{DIM},nvs_data::Ghost_VS_Data,direction::Int,rot::Float64) where{DIM}# Should be merged into gradient_encode?
    heavi = [x<=0. for x in rot.*@views vs_data.midpoint[:,direction]]
    level = vs_data.level;level_n = nvs_data.level
    table_n = zeros(Int,nvs_data.vs_num)
    table_index_n = 1;index = 1;flag = 0.
    for i in 1:vs_data.vs_num
        if level[i]==level_n[index]
            index+=1
            # root_flag += 1/2^(DIM*level[i])
            if !heavi[i]
                table_n[index] = table_index_n
                table_index_n += 1
            end
        elseif level[i]<level_n[index]
            # push!(vg_index,i)
            while flag != 1.0
                # append!(temp_id,search_velocity_gradient_temp(vs_data,i,@views nvs_data.midpoint[index,:],root,nums,global_data))
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                if !heavi[i]
                    table_n[index] = table_index_n
                    table_index_n += 1
                end
                index+=1
            end
            flag = 0.0
        else
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                if !heavi[i] 
                    table_n[index] = table_index_n
                    table_index_n += 1
                end
                index += 1
                flag = 0.0
            end
        end
    end
    index = 1;flag = 0.
    vg_count_n = 1
    heavi_index_n = Int[]
    for i in 1:vs_data.vs_num
        if level[i]==level_n[index]
            index+=1
            # root_flag += 1/2^(DIM*level[i])
        elseif level[i]<level_n[index]
            # push!(vg_index,i)
            while flag != 1.0
                # append!(temp_id,search_velocity_gradient_temp(vs_data,i,@views nvs_data.midpoint[index,:],root,nums,global_data))
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index+=1
            end
            flag = 0.0
        else
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                if !heavi[i] 
                    append!(heavi_index_n,table[nvg[:,vg_count_n]])
                end
                vg_count_n += 1
                index += 1
                flag = 0.0
            end
        end
    end
    heavi_index_n = reshape(heavi_index_n,DIM,:)
    return nothing,heavi_index_n
end
function velocity_gradient_encode!(face::FullFace,global_data::Global_Data)
    #= 
        Actually here_temp is stored in there_data while the there one is stored in here_data. 
        Here and there refer to which velocity space the gradient describes. 
        But where they are used determines where they are stored.
    =#
    here_vg_temp,there_vg_temp = velocity_gradient_encode(face.here_data.vs_data,face.there_data.vs_data,global_data)
    here_vg_heavi, there_vg_heavi = vg_heavi_id_encode(here_vg_temp,there_vg_temp,face.here_data,face.there_data,face.direction,face.rot)
    id = face.rot<0 ? 2*face.direction : 2*face.direction-1; idn = face.rot<0 ? 2*face.direction-1 : 2*face.direction # The index of the neighbor that the gradient describes.
    push!(face.here_data.neighbor.vg[id],VelocityGradient(there_vg_temp,there_vg_heavi))
    push_vg!(here_vg_heavi,here_vg_temp,face.there_data,idn)
end
function velocity_gradient_encode!(face::HangingFace,global_data::Global_Data)
    #= 
        Actually here_temp is stored in there_data while the there one is stored in here_data. 
        Here and there refer to which velocity space the gradient describes. 
        But where they are used determines where they are stored.
    =#
    id = face.rot<0 ? 2*face.direction : 2*face.direction-1; idn = face.rot<0 ? 2*face.direction-1 : 2*face.direction
    for i in eachindex(face.there_data)
        here_vg_temp,there_vg_temp = velocity_gradient_encode(face.here_data.vs_data,face.there_data[i].vs_data,global_data)
        here_vg_heavi, there_vg_heavi = vg_heavi_id_encode(here_vg_temp,there_vg_temp,face.here_data,face.there_data[i],face.direction,face.rot)
        push!(face.here_data.neighbor.vg[id],VelocityGradient(there_vg_temp,there_vg_heavi))
        push_vg!(here_vg_heavi,here_vg_temp,face.there_data[i],idn)
    end
end
function velocity_gradient_encode!(face::BackHangingFace,global_data::Global_Data)
    #= 
        Actually here_temp is stored in there_data while the there one is stored in here_data. 
        Here and there refer to which velocity space the gradient describes. 
        But where they are used determines where they are stored.
    =#
    id = face.rot<0 ? 2*face.direction : 2*face.direction-1; idn = face.rot<0 ? 2*face.direction-1 : 2*face.direction
    for i in eachindex(face.here_data)
        here_vg_temp,there_vg_temp = velocity_gradient_encode(face.here_data[i].vs_data,face.there_data.vs_data,global_data)
        here_vg_heavi, there_vg_heavi = vg_heavi_id_encode(here_vg_temp,there_vg_temp,face.here_data[i],face.there_data,face.direction,face.rot)
        push!(face.here_data[i].neighbor.vg_template[id],VelocityGradient(there_vg_temp,there_vg_heavi))
        push_vg!(here_vg_heavi,here_vg_temp,face.there_data,idn)
    end
end
function search_velocity_gradient_temp(vs_data::AbstractVsData,id,midpoint_n::AbstractVector,root::Array,num::Array,global_data::Global_Data{DIM}) where{DIM}
    # Lookup the temp in a single root.
    midpoints = vs_data.midpoints;midpoint = @views midpoints[id,:];level = vs_data.level[id]
    offset = sign.(midpoint_n-midpoint)
    quadrature = global_data.config.quadrature
    trees_num = global_data.config.vs_trees_num
    ds = [(quadrature[2*i]-quadrature[2*i-1])/trees_num[i] for i in 1:DIM]
    offset_midpoint = midpoint+offset.*ds/2^level
    root_id = [fld(midpoint[i]-quadrature[i],ds[i])+1 for i in 1:DIM]
    temp_root = root_id+offset
    temp_id = Vector{Int}(undef,DIM)
    for i in eachindex(temp_root)
        if temp_root[i]>trees_num[i]||temp_root[i]<1 # Bound check.
            @show temp_root
            throw(`Border velocity cell is refined! A larger velocity domain is required.`)
        end
        if offset_midpoint[i]*midpoint[i]<0 # Zero edge. The search range is required to be local for the upwind flux interpolation.
            search_range = @views midpoints[root[root_id]-num[root_id]+1:root[root_id],:]
            range_offset = root[root_id]-num[root_id]
            _,index = findmin(x->abs(x[i]-midpoint[i])>EPS ? norm(x-midpoint) : Inf,eachrow(search_range)) # Exclude the colinear case along the interested direction.
        else
            if fld(offset_midpoint[i]-quadrature[i],ds[i])+1==root_id[i]
                search_range = @views midpoints[root[root_id]-num[root_id]+1:root[root_id],:]
                range_offset = root[root_id]-num[root_id]
            else
                search_id = copy(root_id);search_id[i] = temp_root[i]
                search_range = @views midpoints[root[search_id]-num[seach_id]+1:root[search_id],:]
                range_offset = root[search_id]-num[search_id]
            end
            _,index = findmin(x->offset[i]*(x[i]-midpoint[i])>EPS ? norm(x-midpoint_n) : Inf,eachrow(search_range))
        end
        temp_id[i] = index+range_offset
    end
    return temp_id
end
function velocity_gradient_encode(vs_data::VS_Data,nvs_data::VS_Data,global_data::Global_Data{DIM}) where{DIM}
    root,nums = velocity_root_encode(vs_data,global_data)
    nroot,nnums = velocity_root_encode(nvs_data,global_data)
    level = vs_data.level;level_n = nvs_data.level
    # vg_index = Int[];vg_index_n = Int[]
    temp_id = Int[];temp_id_n = Int[]
    index = 1;flag = 0.
    # root_flag = 0.; root_flag_n = 0.
    # root_id = root_idn = 1 # The current root indices of two sides for the convenience of the templates search.
    for i in 1:vs_data.vs_num
        if level[i]==level_n[index]
            index+=1
            # root_flag += 1/2^(DIM*level[i])
        elseif level[i]<level_n[index]
            # push!(vg_index,i)
            while flag != 1.0
                @views append!(temp_id,search_velocity_gradient_temp(vs_data,i,nvs_data.midpoint[index,:],root,nums,global_data))
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index+=1
            end
            flag = 0.0
        else
            @views append!(temp_id_n,search_velocity_gradient_temp(nvs_data,index,vs_data.midpoint[i,:],nroot,nnums,global_data))
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
    velocity_gradient_temp = reshape(temp_id,2,:);velocity_gradient_temp_n = reshape(temp_id_n,2,:)
    return velocity_gradient_temp,velocity_gradient_temp_n
end
function velocity_gradient_encode(vs_data::VS_Data,nvs_data::Ghost_VS_Data,global_data::Global_Data{DIM}) where{DIM}
    nroot,nnums = velocity_root_encode(nvs_data,global_data)
    level = vs_data.level;level_n = nvs_data.level
    # vg_index = Int[];vg_index_n = Int[]
    temp_id_n = Int[]
    index = 1;flag = 0.
    # root_flag = 0.; root_flag_n = 0.
    # root_id = root_idn = 1 # The current root indices of two sides for the convenience of the templates search.
    for i in 1:vs_data.vs_num
        if level[i]==level_n[index]
            index+=1
            # root_flag += 1/2^(DIM*level[i])
        elseif level[i]<level_n[index]
            # push!(vg_index,i)
            while flag != 1.0
                # append!(temp_id,search_velocity_gradient_temp(vs_data,i,@views nvs_data.midpoint[index,:],root,nums,global_data))
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index+=1
            end
            flag = 0.0
        else
            @views append!(temp_id_n,search_velocity_gradient_temp(nvs_data,index,vs_data.midpoint[i,:],nroot,nnums,global_data))
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
    velocity_gradient_temp_n = reshape(temp_id_n,2,:)
    return nothing,velocity_gradient_temp_n
end
function velocity_root_encode(vs_data::AbstractVsData,global_data::Global_Data{2}) # Generate a 2-dimensional array whose entries are the tail index of each root. The size of the array is vs_trees_num[1]×vs_trees_num[2]...
    trees_num = global_data.config.vs_trees_num
    root = Vector{Int}(undef,trees_num[1]*trees_num[2])
    DIM = 2
    flag = 0.
    index = 1
    level = vs_data.level
    for i in eachindex(level)
        flag+=1/2^(DIM*(level[i]))
        if flag==1.
            root[index] = i
            flag=0.;index+=1
        end
    end
    nums = [i==1 ? root[i] : root[i]-root[i-1] for i in eachindex(root)]
    return reshape(root,trees_num[1],trees_num[2]), reshape(nums,trees_num[1],trees_num[2]) # i~x, j~y
end
function velocity_root_encode(vs_data::AbstractVsData,global_data::Global_Data{3}) # Generate a 3-dimensional array whose entries are the tail index of each root. The size of the array is vs_trees_num[1]×vs_trees_num[2]...
    trees_num = global_data.config.vs_trees_num
    root = Vector{Int}(undef,trees_num[1]*trees_num[2]*trees_num[3])
    DIM = 3
    flag = 0.
    index = 1
    level = vs_data.level
    for i in eachindex(level)
        flag+=1/2^(DIM*(level[i]))
        if flag==1.
            root[index] = i
            flag=0.;index+=1
        end
    end
    nums = [i==1 ? root[i] : root[i]-root[i-1] for i in eachindex(root)]
    return reshape(root,trees_num[1],trees_num[2],trees_num[3]), reshape(nums,trees_num[1],trees_num[2],trees_num[3]) # i~x, j~y, k~z
end
