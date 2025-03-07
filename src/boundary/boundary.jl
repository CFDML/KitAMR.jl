include("Circle.jl")
include("Vertices.jl")
get_bc(bc::AbstractVector) = bc
function calc_IB_ρw_2D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_2D2F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function calc_IB_ρw_3D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_3D1F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function IB_prim(bc::AbstractVector,::AbstractVector,ρw::Real)
    prim = copy(bc)
    prim[1] = ρw
    return prim
end
function calc_solid_cell_slope!(svdata::AbstractVsData{DIM,NDF},fvdata::VS_Data{DIM,NDF},smid::Vector{Float64},fmid::Vector{Float64},direction::Integer) where{DIM,NDF}
    j = 1
    flag = 0.0
    level = svdata.level
    sdf = @view(svdata.sdf[:,:,direction])
    df = svdata.df
    level_n = fvdata.level
    df_n = fvdata.df
    dx = fmid[direction]-smid[direction]
    @inbounds for i in 1:svdata.vs_num
        if level[i] == level_n[j]
            @views @. sdf[i, :] = df_n[j,:]-df[i,:]
            j += 1
        elseif level[i] < level_n[j]
            @views sdf[i, :] .= -df[i,:]
            while flag != 1.0
                @views @. sdf[i, :] += df_n[j, :]/ 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @views @. sdf[i, :] = df_n[j,:]-df[i,:]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
    sdf./=dx
end
function calc_w0(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,::Global_Data{2,2})
    @views micro_to_macro_2D2F(midpoint[:,1],midpoint[:,2],df[:,1],df[:,2],weight)
end
function calc_qf(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,prim::AbstractVector,::Global_Data{2,2})
    @views heat_flux_2D2F(midpoint[:,1],midpoint[:,2],df[:,1],df[:,2],prim,weight)
end
function calc_boundary_qf(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,fprim::AbstractVector,sprim::AbstractVector,Θ::AbstractVector,::Global_Data{2,2})
    fweight = @. weight*(1.0-Θ)
    sweight = @. weight*Θ
    @views heat_flux_2D2F(midpoint[:,1],midpoint[:,2],df[:,1],df[:,2],fprim,fweight)+heat_flux_2D2F(midpoint[:,1],midpoint[:,2],df[:,1],df[:,2],sprim,sweight)
end
function calc_pressure(midpoint::AbstractMatrix,df::AbstractMatrix,weight::AbstractVector,::Global_Data{2})
    @views pressure_2D(midpoint[:,1],midpoint[:,2],df[:,1],weight)
end
function save_boundary_result!(ib::AbstractBoundary,ps_data,solid_neighbor::SolidNeighbor{DIM,NDF,ID},boundary_results,amr::AMR{DIM,NDF}) where{DIM,NDF,ID}
    global_data = amr.global_data
    vs_data = ps_data.vs_data
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    Θ = heaviside.(vn)
    dir = get_dir(ID)
    ib_point = aux_point[dir]+0.5*(ps_data.midpoint[dir]-solid_cell.midpoint[dir])
    ib_df = @views vs_data.df+vs_data.sdf[:,:,dir]*(ib_point-ps_data.midpoint[dir])
    aux_df = zeros(vs_data.vs_num,NDF)
    vs_interpolate!(ib_df,vs_data.level,ib_point,s_vs_data.df,
        s_vs_data.level,solid_cell.midpoint[dir],aux_df,aux_point[dir],amr)
    ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
    s_prim = IB_prim(ib,aux_point,ρw)
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),s_prim,amr.global_data)
        end
    end
    aux_w = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    aux_prim = get_prim(aux_w,global_data)
    # aux_qf = calc_boundary_qf(vs_data.midpoint,aux_df,vs_data.weight,f_prim,s_prim,Θ,global_data)
    aux_qf = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,aux_prim,global_data)
    aux_p = calc_pressure(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    push!(boundary_results[ps_data.bound_enc].midpoints,aux_point)
    push!(boundary_results[ps_data.bound_enc].ps_solutions,Boundary_PS_Solution(aux_prim,aux_qf,aux_p))
end
function save_boundary_result!(ib::AbstractBoundary,ps_data::PS_Data{DIM,NDF},boundary_results::Vector{Boundary_Solution},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        save_boundary_result!(ib,ps_data,ps_data.neighbor.data[i][1],boundary_results,amr)
    end
end
function solid_cell_index_encoder!(solid_cell_index::Vector{Int},now_index::Int)
    id = findfirst(x->x==0,solid_cell_index)
    isnothing(id) && (@error `A larger SOLID_CELL_ID_NUM is needed!`)
    solid_cell_index[id]=now_index
end
function solid_cell_index_decoder(solid_cell_index::Vector{Int})
    ids = findall(x->x!=0,solid_cell_index)
    return solid_cell_index[ids]
end
function solid_cell_index2ranks(indices::Vector,quadids::Vector,gfq::Vector) # give the 1-based ranks containing solid_cells in indices
    qids = quadids[indices]
    lrank = MPI.Comm_rank(MPI.COMM_WORLD)+1
    ranks = Int[];lids = Int[]
    for i in eachindex(qids)
        rank = findfirst(x->x>qids[i],gfq)-1
        !(rank==lrank||in(rank,ranks))&&push!(ranks,rank)
        rank==lrank&&push!(lids,indices[i])
    end
    return ranks,lids
end
function solid_flag(::Domain{Maxwellian},::AbstractVector) # Does midpoint locate at solid?
    return false
end

function solid_cell_flag(boundaries::Vector{AbstractBoundary},midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data)
    for boundary in boundaries
        solid_cell_flag(boundary,midpoint,ds,global_data,solid_flag(boundary,midpoint))&& return true
    end
    return false
end
function InsideSolid_flag(boundary::AbstractBoundary,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data)
    inside = solid_flag(boundary,midpoint)
    (solid_flag(boundary,midpoint)&&!solid_cell_flag(boundary,midpoint,ds,global_data,inside))&& return true # In solid region and not the ghost cell( a.k.a. solid cell)
    return false
end
function InsideSolid_flag(boundaries::Vector{AbstractBoundary},midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data)
    for boundary in boundaries
        InsideSolid_flag(boundary,midpoint,ds,global_data) && return true
    end
    return false
end

function pre_broadcast_boundary_points(boundary_points::Vector)
    Nb = length(boundary_points)
    numbers = Vector{Int}(undef,Nb)
    for i in eachindex(boundary_points)
        numbers[i] = length(boundary_points[i])
    end
    Numbers = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    for i in eachindex(Numbers)
        if i-1 == MPI.Comm_rank(MPI.COMM_WORLD)
            Numbers[i] = numbers
        else
            Numbers[i] = Vector{Int}(undef,Nb)
        end
    end
    for i in eachindex(Numbers)
        MPI.Bcast!(Numbers[i],i-1,MPI.COMM_WORLD)
    end
    return Numbers
end
function broadcast_boundary_midpoints!(boundary_points::Vector{Vector{Vector{Float64}}},::Global_Data{DIM})where{DIM}
    Numbers = pre_broadcast_boundary_points(boundary_points)
    rbuffer = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points)) # boundaries{ranks{points}}
    sbuffer = Vector{Vector{Float64}}(undef,length(boundary_points)) # boundaries{points}
    for i in eachindex(boundary_points)
        buffer = Vector{Float64}(undef,DIM*length(boundary_points[i]))
        for j in eachindex(boundary_points[i])
            buffer[DIM*(j-1)+1:DIM*j] .= boundary_points[i][j]
        end
        sbuffer[i] = buffer
    end
    for i in eachindex(boundary_points)
        rbuffer[i] = Vector{Vector{Float64}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
        for j in eachindex(Numbers)
            if j-1==MPI.Comm_rank(MPI.COMM_WORLD) 
                rbuffer[i][j] = sbuffer[i]
            else
                rbuffer[i][j] = Vector{Float64}(undef,DIM*Numbers[j][i])
            end
        end
    end
    for i in eachindex(boundary_points)
        for j in eachindex(Numbers)
            MPI.Bcast!(rbuffer[i][j],j-1,MPI.COMM_WORLD)
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
    solid_points_global = Vector{Vector{Vector{Float64}}}(undef,length(boundary_points))
    for i in eachindex(boundary_points)
        solid_points_global[i] = Vector{Float64}[]
        for j in eachindex(Numbers)
            if j-1 == MPI.Comm_rank(MPI.COMM_WORLD)
                append!(solid_points_global[i],boundary_points[i])
            else
                for k in 1:Int(length(rbuffer[i][j])/DIM)
                    push!(solid_points_global[i],rbuffer[i][j][DIM*(k-1)+1:DIM*k]) 
                end
            end
        end
    end
    return solid_points_global
end


function IB_flag(boundary::AbstractBoundary,aux_point::AbstractVector,midpoint::AbstractVector,::AbstractVector)
    r = boundary.search_radius
    norm(midpoint-aux_point)<r
end
function IB_flag(boundaries::Vector{AbstractBoundary},aux_points::Vector{Vector{Vector{Float64}}},midpoint::AbstractVector,ds::AbstractVector)
    for i in eachindex(boundaries)
        solid_flag(boundaries[i],midpoint) && return false
        for j in eachindex(aux_points[i])
            IB_flag(boundaries[i],aux_points[i][j],midpoint,ds) && return true
        end
    end
    return false
end

function vs_extrapolate!(df::AbstractMatrix{Float64},sdf::AbstractMatrix{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},dx::Float64,::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(dft,1)
        if levelt[i] == level[j]
            @. dft[i, :] = @views df[j, :]+sdf[j,:]*dx
            j += 1
        elseif levelt[i] < level[j]
            while flag != 1.0
                @. dft[i, :] += @views (df[j, :]+sdf[j,:]*dx)/ 2^(DIM * (level[j] - levelt[i]))
                flag += 1 / 2^(DIM * (level[j] - levelt[i]))
                j += 1
            end
            flag = 0.0
        else
            @. dft[i, :] = @views df[j,:]+sdf[j,:]*dx
            flag += 1 / 2^(DIM * (levelt[i] - level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function vs_extrapolate!(df::AbstractMatrix{Float64},sdf::AbstractMatrix{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},dx::Float64,weight::AbstractVector{Float64},::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(dft,1)
        if levelt[i] == level[j]
            @. dft[i, :] += @views (df[j, :]+sdf[j,:]*dx)*weight[i]
            j += 1
        elseif levelt[i] < level[j]
            while flag != 1.0
                @. dft[i, :] += @views (df[j, :]+sdf[j,:]*dx)/ 2^(DIM * (level[j] - levelt[i]))*weight[i]
                flag += 1 / 2^(DIM * (level[j] - levelt[i]))
                j += 1
            end
            flag = 0.0
        else
            @. dft[i, :] += @views (df[j,:]+sdf[j,:]*dx)*weight[i]
            flag += 1 / 2^(DIM * (levelt[i] - level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function vs_extrapolate!(df::AbstractMatrix{Float64},sdf::AbstractArray{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},dx::Vector{Float64},weight::AbstractVector{Float64},::AMR{DIM,NDF}) where{DIM,NDF}
    ddf = [dot(sdf[i,j,:],dx) for i in axes(sdf,1), j in axes(sdf,2)]
    j = 1;flag = 0.0
    @inbounds for i in axes(dft,1)
        if levelt[i] == level[j]
            @. dft[i, :] += @views (df[j, :]+ddf[j,:])*weight[i]
            j += 1
        elseif levelt[i] < level[j]
            while flag != 1.0
                @. dft[i, :] += @views (df[j, :]+ddf[j,:])/ 2^(DIM * (level[j] - levelt[i]))*weight[i]
                flag += 1 / 2^(DIM * (level[j] - levelt[i]))
                j += 1
            end
            flag = 0.0
        else
            @. dft[i, :] += @views (df[j,:]+ddf[j,:])*weight[i]
            flag += 1 / 2^(DIM * (levelt[i] - level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function update_solid_cell!(ps_data::PS_Data{2,NDF},fluid_cells::Vector,amr::AMR{2,NDF}) where{NDF}
    vs_data = ps_data.vs_data;vs_data.df.=0.
    weights = Matrix{Float64}(undef,vs_data.vs_num,length(fluid_cells))
    for i in eachindex(fluid_cells)
        l = ps_data.midpoint-fluid_cells[i].midpoint;l/=norm(l)
        weights[:,i] .= [max(0.,dot(u,l)/norm(u))^2 for u in eachrow(vs_data.midpoint)]
    end
    weight_i = Vector{Float64}(undef,vs_data.vs_num)
    weight_sum = sum(weights,dims=2)
    for i in eachindex(fluid_cells)
        f_vs_data = fluid_cells[i].vs_data
        for j in eachindex(weight_i)
            weight_i[j] = weight_sum[j]==0. ? 1.0/length(fluid_cells) : weights[j,i]/weight_sum[j]
        end
        fdf = f_vs_data.df;fsdf = f_vs_data.sdf;dx = ps_data.midpoint-fluid_cells[i].midpoint
        vs_extrapolate!(fdf,fsdf,f_vs_data.level,vs_data.df,vs_data.level,dx,weight_i,amr)
    end
end
function update_solid_cell!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = [ps_data.neighbor.data[i][1] for i in fluid_dirs]
    update_solid_cell!(ps_data,fluid_cells,amr)
end
function update_solid_cell!(amr::AMR)
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc>=0)&&continue
            update_solid_cell!(ps_data,amr)
        end
    end
end
function initialize_solid_neighbor!(amr::AMR)
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            initialize_solid_neighbor!(ps_data,amr)
        end
    end
end

function initialize_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_dirs = findall(x->!isnothing(x[1])&&x[1].bound_enc<0,ps_data.neighbor.data)
    vs_data = ps_data.vs_data
    for i in solid_dirs
        solid_cell = ps_data.neighbor.data[i][1]
        if ps_data.bound_enc!=0
            ps_data.bound_enc!=-solid_cell.bound_enc&&throw(`Different immersed boundaries correspond to the same fluid cell!`)
        end
        ps_data.bound_enc=-solid_cell.bound_enc
        ib = amr.global_data.config.IB[-solid_cell.bound_enc]
        aux_point,normal = calc_intersect(ps_data.midpoint,solid_cell.midpoint,ib)
        ic = get_bc(ib.bc)
        svsdata = VS_Data{DIM,NDF}(
            vs_data.vs_num,
            vs_data.level,
            vs_data.weight,
            vs_data.midpoint,
            discrete_maxwell(vs_data.midpoint,ic,amr.global_data),
            zeros(vs_data.vs_num,NDF,DIM),
            Matrix{Float64}(undef,0,0)
        )
        ps_data.neighbor.data[i][1] = SolidNeighbor{DIM,NDF,i}(
            solid_cell.bound_enc,aux_point,normal,solid_cell,solid_cell.midpoint,svsdata
        )
    end
end
function vs_interpolate!(f_df::AbstractMatrix,f_level::AbstractVector{Int8},fx,s_df,s_level,sx,b_df,bx,::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(f_df,1)
        if f_level[i] == s_level[j]
            @. b_df[i, :] = @views f_df[i,:]+(s_df[j, :]-f_df[i,:])/(sx-fx)*(bx-fx)
            j += 1
        elseif f_level[i] < s_level[j]
            @. b_df[i,:] = @views -f_df[i,:]
            while flag != 1.0
                @. b_df[i, :] += @views s_df[j,:]/ 2^(DIM * (s_level[j] - f_level[i]))
                flag += 1 / 2^(DIM * (s_level[j] - f_level[i]))
                j += 1
            end
            @. b_df[i,:]*=(bx-fx)/(sx-fx)
            @. b_df[i,:] += @views f_df[i,:]
            flag = 0.0
        else
            @. b_df[i, :] += @views f_df[i,:]+(s_df[j, :]-f_df[i,:])/(sx-fx)*(bx-fx)
            flag += 1 / 2^(DIM * (f_level[i] - s_level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF,ID},amr::AMR) where{DIM,NDF,ID}
    ib = amr.global_data.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    dir = get_dir(ID)
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    dxL = 0.5*ps_data.ds[dir]
    dxR = norm(solid_cell.midpoint-aux_point)
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    ib_point = aux_point[dir]+0.5*(ps_data.midpoint[dir]-solid_cell.midpoint[dir])
    ib_df = @views vs_data.df+vs_data.sdf[:,:,dir]*(ib_point-ps_data.midpoint[dir])
    Θ = heaviside.(vn)
    vs_interpolate!(ib_df,vs_data.level,ib_point,s_vs_data.df,
        s_vs_data.level,solid_cell.midpoint[dir],aux_df,aux_point[dir],amr)
    ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
    aux_prim = IB_prim(ib,aux_point,ρw)
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),aux_prim,amr.global_data)
        end
    end
    @. solid_neighbor.vs_data.df = aux_df+(aux_df-ib_df)/dxL*dxR
end
function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        update_solid_neighbor!(ps_data,ps_data.neighbor.data[i][1],amr)
    end
end
function update_solid_neighbor!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
            update_solid_neighbor!(ps_data,amr)
        end
    end
end
function update_solid!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    initialize_solid_neighbor!(amr)
end