include("Circle.jl")
include("Vertices.jl")
include("Period.jl")
include("immersed_boundary.jl")
get_bc(bc::AbstractVector) = copy(bc)
function calc_IB_ρw_2D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_2D2F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function calc_IB_ρw_3D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_3D1F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(midpoint[:,3]),@view(df[:,1]), bc, weight, Θ, vn)
end
function IB_prim(bc::AbstractVector,::AbstractVector,ρw::Real)
    prim = copy(bc)
    prim[1] = ρw
    return prim
end
function calc_solid_cell_slope!(svdata::AbstractVsData{DIM,NDF},fvdata::AbstractVsData{DIM,NDF},smid::Vector{Float64},fmid::Vector{Float64},direction::Integer) where{DIM,NDF}
    j = 1
    flag = 0.0
    level = svdata.level
    sdf = @view(svdata.sdf[:,:,direction])
    # fsdf = @views fvdata.sdf[:,:,direction]
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

function update_ghost_target_cells_slope!(ps_data::Ghost_PS_Data{DIM,NDF},nps_data::PS_Data{DIM,NDF},dir::Int) where{DIM,NDF}
    df = ps_data.vs_data.df; sdf = @views ps_data.vs_data.sdf[:,:,dir]
    ndf = nps_data.vs_data.df;nsdf = @views nps_data.vs_data.sdf[:,:,dir]
    dx = ps_data.midpoint[dir]-nps_data.midpoint[dir]
    level = ps_data.vs_data.level; level_n = nps_data.vs_data.level
    index = 1;flag = 0.
    @inbounds for i in axes(df,1)
        if level[i]==level_n[index]
            for j in axes(df,2)
                sdf[i,j] = 2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j]
            end
            index += 1
        elseif level[i]<level_n[index]
            sdf[i,:] .= 0.
            while flag != 1.0
                for j in axes(df,2)
                    sdf[i,j] += (2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j])/2^(DIM * (level_n[index] - level[i]))
                end
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index += 1
            end 
            flag = 0.
        else
            for j in axes(df,2)
                sdf[i,j] = 2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j]
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
end
# function save_boundary_result!(ib::AbstractBoundary,ps_data,solid_neighbor::SolidNeighbor{DIM,NDF},boundary_results,amr::AMR{DIM,NDF}) where{DIM,NDF}
#     global_data = amr.global_data;ib = global_data.config.IB[ps_data.bound_enc]
#     vs_data = ps_data.vs_data
#     aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
#     faceid = solid_neighbor.faceid;opp = faceid%2==0 ? faceid-1 : faceid+1
#     dir = get_dir(faceid)
#     # rot = get_rot(faceid)
#     # solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
#     # dxL = 0.5*ps_data.ds[dir]
#     # dxR = norm(solid_neighbor.midpoint-aux_point)
#     vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
#     aux_df = zeros(vs_data.vs_num,NDF)
#     # ib_point = aux_point+0.5*(ps_data.midpoint-solid_neighbor.midpoint)
#     # ib_df =  @views vs_data.df+vs_data.sdf[:,:,dir]*(ib_point[dir]-ps_data.midpoint[dir])
#     Θ = heaviside.(vn)
#     # vs_interpolate!(ib_df,vs_data.level,ib_point[dir],solid_neighbor.ex_df,
#         # vs_data.level,solid_neighbor.midpoint[dir],aux_df,aux_point[dir],amr)
#     upwind1st = ps_data.neighbor.data[opp][1]
#     upwind2nd_extrapolate!(ps_data,solid_neighbor,upwind1st,aux_df,Θ,faceid)
#     # ib_df = vs_data.df+(2*aux_point[dir]-solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])*(@views vs_data.sdf[:,:,dir])
#     cvc_gas_correction!(aux_df,solid_neighbor)
#     aux_prim = get_bc(ib.bc);aux_prim[1] = 1.
#     M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
#     Mu_L,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
#     ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
#     aux_prim = IB_prim(ib,aux_point,ρw)
#     F = aux_prim[1]*M
#     for i in 1:vs_data.vs_num
#         if Θ[i]==1.
#             aux_df[i,:] .= @views F[i,:]
#         end
#     end
#     cvc_correction!(aux_df,F,solid_neighbor,amr)

#     #= AP Maxwell Prototype
#     # end
#     # s_prim = IB_prim(ib,aux_point,ρw)
#     # for i in 1:vs_data.vs_num
#     #     if Θ[i]==1.
#     #         aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),s_prim,amr.global_data)
#     #     end
#     # end
#     # cvc_correction!(aux_df,s_prim,vn,solid_neighbor,amr)
#     =#
#     aux_w = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,global_data)
#     aux_prim = get_prim(aux_w,global_data)
#     aux_qf = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,aux_prim,global_data)
#     aux_p = calc_pressure(vs_data.midpoint,aux_df,vs_data.weight,global_data)
#     push!(boundary_results[ps_data.bound_enc].midpoints,aux_point)
#     push!(boundary_results[ps_data.bound_enc].normal,n)
#     push!(boundary_results[ps_data.bound_enc].ps_solutions,Boundary_PS_Solution(aux_prim,aux_qf,aux_p))
#     dir_path = "./boundary_vs"
#     !isdir(dir_path) && mkpath(dir_path)
#     if NDF==2
#         @suppress write_vs_VTK(aux_df,vs_data,amr,dir_path*"/"*string(ps_data.midpoint)*string(n),["h","b"],fieldvalues_fn)
#     elseif NDF==3
#         @suppress write_vs_VTK(aux_df,vs_data,amr,dir_path*"/"*string(ps_data.midpoint)*string(n),["df"],fieldvalues_fn)
#     end
# end
function save_boundary_result!(ib::AbstractBoundary,ps_data,solid_neighbor::SolidNeighbor{DIM,NDF},boundary_results,amr::AMR{DIM,NDF}) where{DIM,NDF}
    global_data = amr.global_data;ib = global_data.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid);rot = get_rot(faceid)
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    opp = faceid%2==0 ? faceid-1 : faceid+1
    upwind1 = ps_data.neighbor.data[opp][1];vs_data1 = upwind1.vs_data
    Θ = heaviside.(vn) # Solid to gas (emission) is 1.
    Θ_grid = @views heaviside.(rot*vs_data.midpoint[:,dir]) # Solid to gas (downwind) is 1.
    construct_incident_df!(aux_df,ps_data,solid_neighbor;Θ,Θ_grid)
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
    Mu_L,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
    ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
    aux_prim = IB_prim(ib,aux_point,ρw)
    F = aux_prim[1]*M
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            aux_df[i,:] .= F[i,:]
        end
    end
    cvc_correction!(aux_df,F,solid_neighbor,amr)
    aux_w = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    aux_prim = get_prim(aux_w,global_data)
    aux_qf = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,aux_prim,global_data)
    aux_p = calc_pressure(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    push!(boundary_results[ps_data.bound_enc].midpoints,aux_point)
    push!(boundary_results[ps_data.bound_enc].normal,n)
    push!(boundary_results[ps_data.bound_enc].ps_solutions,Boundary_PS_Solution(aux_prim,aux_qf,aux_p))
    dir_path = "./boundary_vs"
    !isdir(dir_path) && mkpath(dir_path)
    if NDF==2
        @suppress write_vs_VTK(aux_df,vs_data,amr,dir_path*"/"*string(ps_data.midpoint)*string(n),["h","b"],fieldvalues_fn)
    elseif NDF==3
        @suppress write_vs_VTK(aux_df,vs_data,amr,dir_path*"/"*string(ps_data.midpoint)*string(n),["df"],fieldvalues_fn)
    end
end
function save_boundary_result!(ib::AbstractBoundary,ps_data::PS_Data{DIM,NDF},boundary_results::Vector{Boundary_Solution},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        save_boundary_result!(ib,ps_data,ps_data.neighbor.data[i][1],boundary_results,amr)
    end
end
# function solid_cell_index_encoder!(solid_cell_index::Vector{Int},now_index::Int)
#     id = findfirst(x->x==0,solid_cell_index)
#     isnothing(id) && (@error `A larger SOLID_CELL_ID_NUM is needed!`)
#     solid_cell_index[id]=now_index
# end
# function solid_cell_index_decoder(solid_cell_index::Vector{Int})
#     ids = findall(x->x!=0,solid_cell_index)
#     return solid_cell_index[ids]
# end
# function solid_cell_index2ranks(indices::Vector,quadids::Vector,gfq::Vector) # give the 1-based ranks containing solid_cells in indices
#     qids = quadids[indices]
#     lrank = MPI.Comm_rank(MPI.COMM_WORLD)+1
#     ranks = Int[];lids = Int[]
#     for i in eachindex(qids)
#         rank = findfirst(x->x>qids[i],gfq)-1
#         !(rank==lrank||in(rank,ranks))&&push!(ranks,rank)
#         rank==lrank&&push!(lids,indices[i])
#     end
#     return ranks,lids
# end
function solid_flag(::Domain,::AbstractVector) # Does midpoint locate at solid?
    return false
end
function solid_flag(midpoint,global_data::Global_Data)
    for boundary in global_data.config.IB
        solid_flag(boundary,midpoint) && return true
    end
    for domain in global_data.config.domain
        solid_flag(domain,midpoint) && return true
    end
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

# function pre_broadcast_boundary_points(boundary_points::Vector)
#     Nb = length(boundary_points)
#     numbers = Vector{Int}(undef,Nb)
#     for i in eachindex(boundary_points)
#         numbers[i] = length(boundary_points[i])
#     end
#     Numbers = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
#     for i in eachindex(Numbers)
#         if i-1 == MPI.Comm_rank(MPI.COMM_WORLD)
#             Numbers[i] = numbers
#         else
#             Numbers[i] = Vector{Int}(undef,Nb)
#         end
#     end
#     for i in eachindex(Numbers)
#         MPI.Bcast!(Numbers[i],i-1,MPI.COMM_WORLD)
#     end
#     return Numbers
# end

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

function IB_flag(boundary::AbstractBoundary,aux_point::AbstractVector,midpoint::AbstractVector,ds::AbstractVector,level::Integer,global_data::Global_Data{DIM}) where{DIM}
    r = boundary.search_radius
    if level==global_data.config.solver.AMR_PS_MAXLEVEL
        norm(midpoint-aux_point)<r+EPS
    else
        @inbounds  for i in 1:2^DIM
            norm(midpoint+ds.*RMT[DIM][i]-aux_point)<r+EPS && return true
        end
        return false
    end
end
function IB_flag(aux_points::Vector{Vector{Vector{Float64}}},midpoint::AbstractVector,ds::AbstractVector,level::Int8,global_data)
    boundaries = global_data.config.IB
    for i in eachindex(boundaries)
        solid_flag(boundaries[i],midpoint) && return false
        for j in eachindex(aux_points[i])
            IB_flag(boundaries[i],aux_points[i][j],midpoint,ds,level,global_data) && return true
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
    ddf = @views [dot(sdf[i,j,:],dx) for i in axes(sdf,1), j in axes(sdf,2)]
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
# function update_solid_cell!(::Union{DVM,CAIDVM,UGKS},ps_data::PS_Data{DIM,NDF},fluid_cells::Vector,amr::AMR{DIM,NDF}) where{DIM,NDF}
#     vs_data = ps_data.vs_data;vs_data.df.=0.
#     weights = Matrix{Float64}(undef,vs_data.vs_num,length(fluid_cells))
#     for i in eachindex(fluid_cells)
#         l = ps_data.midpoint-fluid_cells[i].midpoint;l/=norm(l)
#         weights[:,i] .= [max(0.,dot(u,l)/norm(u))^2 for u in eachrow(vs_data.midpoint)]
#     end
#     weight_i = Vector{Float64}(undef,vs_data.vs_num)
#     weight_sum = sum(weights,dims=2)
#     for i in eachindex(fluid_cells)
#         f_vs_data = fluid_cells[i].vs_data
#         for j in eachindex(weight_i)
#             weight_i[j] = weight_sum[j]==0. ? 1.0/length(fluid_cells) : weights[j,i]/weight_sum[j]
#         end
#         fdf = f_vs_data.df;fsdf = f_vs_data.sdf;dx = ps_data.midpoint-fluid_cells[i].midpoint
#         vs_extrapolate!(fdf,fsdf,f_vs_data.level,vs_data.df,vs_data.level,dx,weight_i,amr)
#     end
#     ps_data.w = calc_w0(ps_data)
#     ps_data.prim = get_prim(ps_data,amr.global_data)
# end
function update_ex_df!(::Union{DVM,CAIDVM,UGKS},ps_data::PS_Data{DIM,NDF},sn::SolidNeighbor{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    vs_data = ps_data.vs_data
    corner_neighbors = [x[1] for x in ps_data.neighbor.data]
    fluid_dirs = findall(x->!(isa(x,AbstractInsideSolidData)||x.bound_enc<0)&&dot(x.midpoint-ps_data.midpoint,sn.midpoint-ps_data.midpoint)>-EPS,corner_neighbors)
    # fluid_dirs = findall(x->!(isa(x,AbstractInsideSolidData)||x.bound_enc<0),corner_neighbors)
    sn.ex_df.=0.;weights = Matrix{Float64}(undef,vs_data.vs_num,length(fluid_dirs)+1)
    for i in eachindex(fluid_dirs)
        fluid_cell = corner_neighbors[fluid_dirs[i]]
        l = ps_data.midpoint-fluid_cell.midpoint;l/=norm(l)
        weights[:,i] .= [max(0.,dot(u,l)/norm(u))^2 for u in eachrow(vs_data.midpoint)]
    end
    l = sn.midpoint-ps_data.midpoint;l/=norm(l)
    weights[:,end] .= [max(0.,dot(u,l)/norm(u))^2 for u in eachrow(vs_data.midpoint)]


    weight_i = Vector{Float64}(undef,vs_data.vs_num)
    weight_sum = sum(weights,dims=2)
    for i in eachindex(fluid_dirs)
        fluid_cell = corner_neighbors[fluid_dirs[i]]
        f_vs_data = fluid_cell.vs_data
        for j in eachindex(weight_i)
            weight_i[j] = weight_sum[j]==0. ? 1.0/(length(fluid_dirs)+1) : weights[j,i]/weight_sum[j]
        end
        fdf = f_vs_data.df;fsdf = f_vs_data.sdf;dx = sn.midpoint-fluid_cell.midpoint
        # if amr.global_data.status.sim_time>4.0
        #     error_index = findall(x->x==0.,fsdf);error_num = length(error_index)
        #     @show error_num
        # end
        vs_extrapolate!(fdf,fsdf,f_vs_data.level,sn.ex_df,vs_data.level,dx,weight_i,amr)
    end
    f_vs_data = ps_data.vs_data
    for j in eachindex(weight_i)
        weight_i[j] = weight_sum[j]==0. ? 1.0/(length(fluid_dirs)+1) : weights[j,end]/weight_sum[j]
    end
    fdf = f_vs_data.df;fsdf = f_vs_data.sdf;dx = sn.midpoint-ps_data.midpoint
    vs_extrapolate!(fdf,fsdf,f_vs_data.level,sn.ex_df,vs_data.level,dx,weight_i,amr)
end
function update_ex_df!(ps_data,amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_dirs = findall(x->!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc<0,ps_data.neighbor.data[1:2*DIM])
    for dir in solid_dirs
        update_ex_df!(amr.global_data.config.solver.flux,ps_data,ps_data.neighbor.data[dir][1],amr)
    end
end
function update_ex_df!(amr::AMR)
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
            update_ex_df!(ps_data,amr)
        end
    end
end
# function update_solid_cell!(::UGKS,ps_data::PS_Data{2,NDF},fluid_cells::Vector,amr::AMR{2,NDF}) where{NDF}
#     vs_data = ps_data.vs_data;vs_data.df.=0.;ps_data.w .= 0.
#     weights = Matrix{Float64}(undef,vs_data.vs_num,length(fluid_cells))
#     for i in eachindex(fluid_cells)
#         l = ps_data.midpoint-fluid_cells[i].midpoint;l/=norm(l)
#         weights[:,i] .= [max(0.,dot(u,l)/norm(u))^2 for u in eachrow(vs_data.midpoint)]
#     end
#     weight_i = Vector{Float64}(undef,vs_data.vs_num)
#     weight_sum = sum(weights,dims=2)
#     for i in eachindex(fluid_cells)
#         f_vs_data = fluid_cells[i].vs_data
#         for j in eachindex(weight_i)
#             weight_i[j] = weight_sum[j]==0. ? 1.0/length(fluid_cells) : weights[j,i]/weight_sum[j]
#         end
#         fdf = f_vs_data.df;fsdf = f_vs_data.sdf;dx = ps_data.midpoint-fluid_cells[i].midpoint
#         vs_extrapolate!(fdf,fsdf,f_vs_data.level,vs_data.df,vs_data.level,dx,weight_i,amr)
#         dw = @views [dot(sw[i,:],dx) for i in axes(sw,1)]
#         ps_data.w.+=(fluid_cells[i].w+dw)*weight_i
#     end
# end
# function update_solid_cell!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
#     fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
#     fluid_cells = [ps_data.neighbor.data[i][1] for i in fluid_dirs]
#     update_solid_cell!(amr.global_data.config.solver.flux,ps_data,fluid_cells,amr)
# end
# function update_solid_cell!(amr::AMR)
#     for tree in amr.field.trees.data
#         for ps_data in tree
#             (isa(ps_data,InsideSolidData)||ps_data.bound_enc>=0)&&continue
#             update_solid_cell!(ps_data,amr)
#         end
#     end
# end
# function initialize_solid_neighbor!(amr::AMR)
#     for tree in amr.field.trees.data
#         for ps_data in tree
#             (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
#             initialize_solid_neighbor!(ps_data,amr)
#         end
#     end
# end
function initialize_cutted_velocity_cell(n::Vector{Float64},vs_data::VS_Data{2},amr::AMR{2,NDF}) where{NDF}
    any(x->abs(x)<1e-6,n)&&(return CuttedVelocityCells(Int[],Float64[],Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Float64[],Float64[]))
    global_data = amr.global_data
    du = [(global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
        global_data.config.vs_trees_num[i] for i in 1:2]
    vertices = [zeros(2) for _ in 1:4]
    index = Int[];solid_weights = Float64[];gas_weights = Float64[]
    for i in 1:vs_data.vs_num
        ddu = du./2^(vs_data.level[i])
        @views abs(dot(vs_data.midpoint[i,:],n))>0.5*norm(ddu)&&continue
        for j in eachindex(vertices)
            vertices[j] .= @views 0.5*ANTIVT[2][j].*ddu+vs_data.midpoint[i,:]
        end
        flag,gas_weight,solid_weight = cut_rect(n,vertices)
        if flag
            push!(index,i);push!(solid_weights,solid_weight);push!(gas_weights,gas_weight)
        end
    end
    N = length(index)
    gas_dfs = zeros(N,NDF);solid_dfs = zeros(N,NDF)
    weight = copy(vs_data.weight)
    weight[index].=0.
    return CuttedVelocityCells(index,weight,gas_dfs,solid_dfs,gas_weights,solid_weights)
end
function initialize_cutted_velocity_cell(n::Vector{Float64},vs_data::VS_Data{3},amr::AMR{3,NDF}) where{NDF}
    l = findall(x->abs(x)<1e-6,n)
    length(l)>1 &&
        (return CuttedVelocityCells(Int[],Float64[],Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Float64[],Float64[]))
    global_data = amr.global_data
    du = [(global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
        global_data.config.vs_trees_num[i] for i in 1:3]
    vertices = zeros(3,8)
    index = Int[];solid_weights = Float64[];gas_weights = Float64[]
    C,A = cut_cube_CA(n)
    for i in 1:vs_data.vs_num
        ddu = du./2^(vs_data.level[i])
        midpoint = vs_data.midpoint[i,:]
        abs(dot(midpoint,n))>0.5*norm(ddu)&&continue
        for j in axes(vertices,2)
            vertices[:,j] .= 0.5*ANTIVT[3][j].*ddu+midpoint
        end
        flag,gas_weight,solid_weight = cut_cube(n,C,A,midpoint,vertices)
        if flag
            push!(index,i);push!(solid_weights,solid_weight);push!(gas_weights,gas_weight)
        end
    end
    N = length(index)
    gas_dfs = zeros(N,NDF);solid_dfs = zeros(N,NDF)
    weight = copy(vs_data.weight)
    weight[index].=0.
    return CuttedVelocityCells(index,weight,gas_dfs,solid_dfs,gas_weights,solid_weights)
end
function reinit_ib_vs!(amr::AMR)
    trees = amr.field.trees.data
    for i in eachindex(trees)
        for j in eachindex(trees[i])
            ps_data = trees[i][j]
            isa(ps_data,InsideSolidData)&&continue
            if ps_data.bound_enc>0
                ib = amr.global_data.config.IB[ps_data.bound_enc]
                ps_data.prim = get_bc(ib.bc)
                ps_data.w = get_conserved(ps_data,amr.global_data)
                ps_data.vs_data.df = discrete_maxwell(ps_data,amr.global_data)
            end
        end
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
function vs_implicit_interpolate!(f_df::AbstractMatrix,f_level::AbstractVector{Int8},fx,s_df,s_level,sx,b_df,bx,vn,amr::AMR{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0;Δt=amr.global_data.status.Δt
    dx = (bx-fx).-vn*Δt
    @inbounds for i in axes(f_df,1)
        if f_level[i] == s_level[j]
            @. b_df[i, :] = @views f_df[i,:]+(s_df[j, :]-f_df[i,:])/(sx-fx)*dx[i]
            j += 1
        elseif f_level[i] < s_level[j]
            @. b_df[i,:] = @views -f_df[i,:]
            while flag != 1.0
                @. b_df[i, :] += @views s_df[j,:]/ 2^(DIM * (s_level[j] - f_level[i]))
                flag += 1 / 2^(DIM * (s_level[j] - f_level[i]))
                j += 1
            end
            @. b_df[i,:]*=dx[i]/(sx-fx)
            @. b_df[i,:] += @views f_df[i,:]
            flag = 0.0
        else
            @. b_df[i, :] += @views f_df[i,:]+(s_df[j, :]-f_df[i,:])/(sx-fx)*dx[i]
            flag += 1 / 2^(DIM * (f_level[i] - s_level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function cvc_gas_correction!(aux_df,solid_neighbor::SolidNeighbor{DIM,NDF}) where{DIM,NDF}
    cvc = solid_neighbor.cvc
    for i in eachindex(cvc.indices)
        cvc.gas_dfs[i,:] .= @views aux_df[cvc.indices[i],:]
    end
end
function cvc_density(aux_df,ib::AbstractBoundary,vn,Θ,solid_neighbor::SolidNeighbor{2,2})
    vs_data = solid_neighbor.vs_data
    cvc = solid_neighbor.cvc
    n = solid_neighbor.normal
    @inbounds @views SF = sum(@. cvc.weight * vn * aux_df[:,1] * (1.0 - Θ))
    prim = get_bc(ib.bc)
    @inbounds SG = prim[4] / π *
        sum(@. @views cvc.weight * vn * exp(-prim[4] * ((vs_data.midpoint[:,1] - prim[2])^2 + (vs_data.midpoint[:,2] - prim[3])^2)) * Θ)
    for i in eachindex(cvc.indices)
        cvc_vn = vn[cvc.indices[i]]
        SF += cvc.gas_weights[i]*cvc_vn*cvc.gas_dfs[i,1]
        solid_mid = @views vs_data.midpoint[cvc.indices[i],:]
        SG += @views prim[4] / π *cvc.solid_weights[i]*cvc_vn*exp(-prim[4] * ((solid_mid[1] - prim[2])^2 + (solid_mid[2] - prim[3])^2))
    end
    return -SF/SG
end
function cvc_density(aux_df,vn::AbstractVector,Θ,solid_neighbor::SolidNeighbor,MuR::Real)
    cvc = solid_neighbor.cvc
    @inbounds @views SF = sum(@. cvc.weight * vn * aux_df[:,1] * (1.0 - Θ))
    for i in eachindex(cvc.indices)
        cvc_vn = vn[cvc.indices[i]]
        SF += cvc.gas_weights[i]*cvc_vn*cvc.gas_dfs[i,1]
    end
    return -SF/MuR
end
# function cvc_Mu(prim::AbstractVector,vn,Θ,solid_neighbor)
#     vs_data = solid_neighbor.vs_data
#     cvc = solid_neighbor.cvc
#     @inbounds M = @. @views cvc.weight * vn * exp(-prim[4] * ((vs_data.midpoint[:,1] - prim[2])^2 + (vs_data.midpoint[:,2] - prim[3])^2)) 
#     MuL = prim[4]/π*sum(@. M*(1.0-Θ))
#     MuR = prim[4] / π *sum(@. M * Θ)
#     for i in eachindex(cvc.indices)
#         cvc_vn = vn[cvc.indices[i]]
#         solid_mid = @views vs_data.midpoint[cvc.indices[i],:]
#         MuL += @views prim[4] / π *cvc.gas_weights[i]*cvc_vn*exp(-prim[4] * ((solid_mid[1] - prim[2])^2 + (solid_mid[2] - prim[3])^2))
#         MuR += @views prim[4] / π *cvc.solid_weights[i]*cvc_vn*exp(-prim[4] * ((solid_mid[1] - prim[2])^2 + (solid_mid[2] - prim[3])^2))
#     end
#     return MuL,MuR
# end
function cvc_Mu(M::AbstractVector,vn,Θ,solid_neighbor)
    cvc = solid_neighbor.cvc
    Mu_cvc = @. @views cvc.weight * vn *M
    MuL = sum(@. Mu_cvc*(1.0-Θ))
    MuR = sum(@. Mu_cvc * Θ)
    for i in eachindex(cvc.indices)
        cvc_vn = vn[cvc.indices[i]]
        MuL += cvc.gas_weights[i]*cvc_vn*M[cvc.indices[i]]
        MuR += cvc.solid_weights[i]*cvc_vn*M[cvc.indices[i]]
    end
    return MuL,MuR
end

function cvc_correction!(aux_df,aux_prim,solid_neighbor,amr)
    cvc = solid_neighbor.cvc
    vs_data = solid_neighbor.vs_data
    for i in eachindex(cvc.indices)
        cvc.solid_dfs[i,:] .= discrete_maxwell(@view(vs_data.midpoint[cvc.indices[i],:]),aux_prim,amr.global_data)
        @views @. aux_df[cvc.indices[i],:] = (cvc.gas_weights[i]*cvc.gas_dfs[i,:]+cvc.solid_weights[i]*cvc.solid_dfs[i,:])/(cvc.gas_weights[i]+cvc.solid_weights[i])
    end
end
function cvc_correction!(aux_df,F::AbstractMatrix,A::Real,solid_neighbor,amr)
    cvc = solid_neighbor.cvc
    for i in eachindex(cvc.indices)
        cvc.solid_dfs[i,:] .= (1.0+0.5*amr.global_data.status.Δt*A).*@views F[cvc.indices[i],:]
        @views @. aux_df[cvc.indices[i],:] = (cvc.gas_weights[i]*cvc.gas_dfs[i,:]+cvc.solid_weights[i]*cvc.solid_dfs[i,:])/(cvc.gas_weights[i]+cvc.solid_weights[i])
    end
end
function cvc_correction!(aux_df,F::AbstractMatrix,solid_neighbor,amr)
    cvc = solid_neighbor.cvc
    for i in eachindex(cvc.indices)
        cvc.solid_dfs[i,:] .= @views F[cvc.indices[i],:]
        @views @. aux_df[cvc.indices[i],:] = (cvc.gas_weights[i]*cvc.gas_dfs[i,:]+cvc.solid_weights[i]*cvc.solid_dfs[i,:])/(cvc.gas_weights[i]+cvc.solid_weights[i])
    end
end
# function update_solid_neighbor!(::AbstractFluxType,ps_data::PS_Data{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF},amr::AMR) where{DIM,NDF}
#     global_data = amr.global_data;ib = global_data.config.IB[ps_data.bound_enc]
#     vs_data = ps_data.vs_data
#     aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
#     faceid = solid_neighbor.faceid
#     dir = get_dir(faceid)
#     # rot = get_rot(faceid)
#     # solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
#     dxL = 0.5*ps_data.ds[dir]
#     dxR = norm(solid_neighbor.midpoint-aux_point)
#     vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
#     aux_df = zeros(vs_data.vs_num,NDF)
#     ib_point = aux_point+0.5*(ps_data.midpoint-solid_neighbor.midpoint)
#     ib_df =  @views vs_data.df+vs_data.sdf[:,:,dir]*(ib_point[dir]-ps_data.midpoint[dir])
#     Θ = heaviside.(vn)
#     vs_interpolate!(ib_df,vs_data.level,ib_point[dir],solid_neighbor.ex_df,
#         vs_data.level,solid_neighbor.midpoint[dir],aux_df,aux_point[dir],amr)
    
#     # vs_implicit_interpolate!(ib_df,vs_data.level,ib_point[dir],s_vs_data.df,
#     #     s_vs_data.level,solid_cell.midpoint[dir],aux_df,aux_point[dir],vn,amr)
#     # stability_correction!(aux_df,ps_data,ps_data.neighbor.data[Int(faceid+off)][1],dxL,dxR,dir,global_data)
#     # for i in eachindex(Θ)
#     #     if Θ==0&&vs_data.midpoint[i,dir]*rot>0
#     #         @views aux_df[i,:].=vs_data.df[i,:]
#     #         # vs_data.sdf[i,:,dir].= 0.
#     #     end
#     # end
#     cvc_gas_correction!(aux_df,solid_neighbor)
#     aux_prim = get_bc(ib.bc);aux_prim[1] = 1.
#     M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
#     Mu_L,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
#     ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
#     # ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
#     aux_prim = IB_prim(ib,aux_point,ρw)
#     # w0 = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,global_data)
#     # prim0 = get_prim(w0,global_data)
#     # for i in 1:vs_data.vs_num
#     #     if Θ[i]==1.
#     #         aux_df[i,:] .= discrete_maxwell(@view(vs_data.midpoint[i,:]),aux_prim,amr.global_data)
#     #     end
#     # end
#     # cvc_correction!(aux_df,aux_prim,vn,solid_neighbor,amr)
#     F = aux_prim[1]*M
#     # if 1/prim0[end]>1e-6
#     #     gas = global_data.config.gas
#     #     τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
#     #     qf0 = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,prim0,global_data)
#     #     F⁺ = shakhov_part(vs_data.midpoint,F,aux_prim,qf0,global_data)
#     #     Mt = time_int(τ0, global_data.status.Δt)
#     #     sdf = similar(vs_data.sdf)
#     #     sdf[:,:,dir] .= (s_vs_data.df-vs_data.df)/(solid_cell.midpoint[dir]-ps_data.midpoint[dir])
#     #     sdf[:,:,dir%DIM+1] .= @views vs_data.sdf[:,:,dir%DIM+1]
#     #     ddf = [@views dot(sdf[i,j,:],vs_data.midpoint[i,:]) for i in axes(sdf,1),j in axes(sdf,2)]
#     #     # rhs = (global_data.status.Δt-Mt[1])/√(π*aux_prim[end])*0.5*aux_prim[1]+sum(vs_data.weight.*(Mt[1]*@view(F⁺[:,1])+Mt[4]*@view(aux_df[:,1])-Mt[5]*@view(ddf[:,1])).*vn.*(1.0.-Θ))
#     #     rhs = aux_prim[1]*(global_data.status.Δt*Mu_R+Mt[1]*Mu_L)+sum(vs_data.weight.*(Mt[1]*@view(F⁺[:,1])+Mt[4]*@view(aux_df[:,1])-Mt[5]*@view(ddf[:,1])).*vn.*(1.0.-Θ))
#     #     lhs = aux_prim[1]*(-0.5*global_data.status.Δt^2*Mu_R-Mt[3]*Mu_L)
#     #     # lhs = (-0.5*global_data.status.Δt^2+Mt[3])*0.5*aux_prim[1]/√(π*aux_prim[end])
#     #     A = rhs/lhs
#     #     aux_df = (Mt[1]*(F+F⁺)+Mt[3]*A*F+Mt[4]*aux_df-Mt[5]*ddf)/global_data.status.Δt
#     #     cvc_gas_correction!(aux_df,solid_neighbor)
#     #     for i in 1:vs_data.vs_num
#     #         if Θ[i]==1.
#     #             aux_df[i,:] .= (1+0.5*global_data.status.Δt*A).*F[i,:]
#     #         end
#     #     end
#     #     cvc_correction!(aux_df,F,A,solid_neighbor,amr)
#     # else
#         for i in 1:vs_data.vs_num
#             if Θ[i]==1.
#                 aux_df[i,:] .= F[i,:]
#             end
#         end
#         cvc_correction!(aux_df,F,solid_neighbor,amr)
#     # end
#     @. solid_neighbor.vs_data.df = aux_df+(aux_df-ib_df)/dxL*dxR
#     solid_neighbor.w = calc_w0(vs_data.midpoint,solid_neighbor.vs_data.df,vs_data.weight,global_data)
#     solid_neighbor.sw[:,dir] .= (solid_neighbor.w-ps_data.w)./(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])
# end

function update_target_cells_slope!(ps_data::PS_Data{DIM,NDF},nps_data::AbstractPsData,Θ::Vector{Bool},faceid::Int) where{DIM,NDF}
    dir = get_dir(faceid)
    df = ps_data.vs_data.df; sdf = @views ps_data.vs_data.sdf[:,:,dir]
    ndf = nps_data.vs_data.df;nsdf = @views nps_data.vs_data.sdf[:,:,dir]
    dx = ps_data.midpoint[dir]-nps_data.midpoint[dir]
    level = ps_data.vs_data.level; level_n = nps_data.vs_data.level
    index = 1;flag = 0.
    @inbounds for i in axes(df,1)
        if level[i]==level_n[index]
            for j in axes(df,2)
                # s1 = (df[i,j]-ndf[index,j])/dx
                # s2 = s1-nsdf[index,j]
                # s = 2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j]
                # s = s1+s2
                # sdf[i,j] = s1+0.5*(sign(s)+sign(nsdf[index,j]))*sign(s)*s2
                sdf[i,j] = Θ[i] ? 2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j] : (df[i,j]-ndf[index,j])/dx
            end
            index += 1
        elseif level[i]<level_n[index]
            sdf[i,:] .= 0.
            while flag != 1.0
                for j in axes(df,2)
                    sdf[i,j] += (2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j])/2^(DIM * (level_n[index] - level[i]))
                end
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index += 1
            end 
            flag = 0.
        else
            for j in axes(df,2)
                sdf[i,j] = 2.0*(df[i,j]-ndf[index,j])/dx-nsdf[index,j]
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
end
# function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
#     solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
#     for i in solid_neighbors
#         # opp = i%2==0 ? i-1 : i+1
#         # update_target_cells_slope!(ps_data,ps_data.neighbor.data[opp][1],i)
#         update_solid_neighbor!(amr.global_data.config.solver.flux,ps_data,ps_data.neighbor.data[i][1],amr)
#     end
# end
# function update_solid_neighbor!(amr::AMR{DIM,NDF};buffer_steps::Int=0,i::Int = typemax(Int64)) where{DIM,NDF}
#     for tree in amr.field.trees.data
#         for ps_data in tree
#             (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
#             update_solid_neighbor!(ps_data,amr)
#         end
#     end
# end
function update_solid!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    initialize_solid_neighbor!(amr)
    initialize_upwind2nd!(amr)
end
# function update_solid!(amr::AMR{DIM,NDF}) where{DIM,NDF}
#     reinit_solid_neighbor!(amr)
# end
# function initialize_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
#     solid_dirs = findall(x->!isnothing(x[1])&&x[1].bound_enc<0,ps_data.neighbor.data[1:2*DIM])
#     vs_data = ps_data.vs_data
#     for i in solid_dirs
#         solid_cell = ps_data.neighbor.data[i][1]
#         if ps_data.bound_enc!=0
#             ps_data.bound_enc!=-solid_cell.bound_enc&&throw(`Different immersed boundaries correspond to the same fluid cell!`)
#         end
#         ps_data.bound_enc=-solid_cell.bound_enc
#         ib = amr.global_data.config.IB[-solid_cell.bound_enc]
#         aux_point,normal = calc_intersect(ps_data.midpoint,solid_cell.midpoint,ib)
#         svsdata = VS_Data{DIM,NDF}(
#             vs_data.vs_num,
#             vs_data.level,
#             vs_data.weight,
#             vs_data.midpoint,
#             zeros(vs_data.vs_num,NDF),
#             zeros(vs_data.vs_num,NDF,DIM),
#             Matrix{Float64}(undef,0,0)
#         )
#         cvc = initialize_cutted_velocity_cell(normal,svsdata,amr) # heavy overhead
#         ps_data.neighbor.data[i][1] = SolidNeighbor{DIM,NDF}(
#             solid_cell.bound_enc,i,aux_point,normal,solid_cell,solid_cell.midpoint,solid_cell.ds,zeros(DIM+2),zeros(DIM+2,DIM),zeros(vs_data.vs_num,NDF),cvc,svsdata
#         )
#     end
#     return nothing
# end
function initialize_upwind2nd!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    put = Pre_Upwind2nd_Transport()
    ut = Upwind2nd_Transport()
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            if ps_data.bound_enc>0
                solid_dirs = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
                ghost_dirs = findall(x->isa(x[1],Ghost_PS_Data),ps_data.neighbor.data[1:2*DIM])
                for i in ghost_dirs
                    opp = i%2==0 ? i-1 : i+1
                    flag = isa(ps_data.neighbor.data[opp][1],SolidNeighbor)
                    push!(put.mirror_datas[ps_data.neighbor.data[i][1].owner_rank+1],flag)
                end
                for i in solid_dirs
                    opp = i%2==0 ? i-1 : i+1
                    upwind1st = ps_data.neighbor.data[opp][1]
                    solid_neighbor = ps_data.neighbor.data[i][1]
                    if isa(upwind1st,AbstractGhostPsData)
                        push!(ut.solid_neighbors,solid_neighbor) # The solid_neighbors requiring the upwind2nd df.
                        push!(ut.ghost_ps_datas,upwind1st)
                        push!(ut.ranks_sn_ids[upwind1st.owner_rank+1],length(ut.solid_neighbors))
                        push!(ut.ghost_ids[upwind1st.owner_rank+1],upwind1st.ghost_id*DIM+get_dir(i))
                    else
                        upwind2nd = upwind1st.neighbor.data[opp][1]
                        df2 = upwind2nd.vs_data.df
                        level2 = upwind2nd.vs_data.level
                        vs_project!(df2,level2,solid_neighbor.upwind2nd_df,solid_neighbor.vs_data.level,upwind2nd.vs_data)
                    end
                end
            end
            target_dirs = findall(x->!isa(x[1],Nothing)&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>0,ps_data.neighbor.data[1:2*DIM])
            for i in target_dirs
                opp = i%2==0 ? i-1 : i+1
                target_cell = ps_data.neighbor.data[i][1]
                upwind2nd = ps_data.neighbor.data[opp][1]
                if isa(target_cell,Ghost_PS_Data)
                    push!(ut.mirror_ps_datas,upwind2nd) # The order of the mirror data is consistent with upwind1st's ghost_id.
                    push!(ut.target_ps_datas,target_cell)
                    push!(put.ghost_ids[target_cell.owner_rank+1],target_cell.ghost_id*DIM+get_dir(i))
                    push!(put.ranks_ghost_ids[target_cell.owner_rank+1],length(ut.mirror_ps_datas))
                end
            end
            # end
        end
    end
    for i in 1:MPI.Comm_size(MPI.COMM_WORLD)
        if isempty(put.ranks_ghost_ids[i])
            put.ghost_datas[i] = Bool[]
        else
            put.ghost_datas[i] = Vector{Bool}(undef,length(put.ranks_ghost_ids[i]))
        end
    end
    data_exchange!(put.mirror_datas,put.ghost_datas)
    for i in 1:MPI.Comm_size(MPI.COMM_WORLD)
        isempty(put.ranks_ghost_ids[i])&&continue
        ghost_ids = sortperm(put.ghost_ids[i])
        for j in eachindex(ghost_ids)
            if !put.ghost_datas[i][j]
                put.ranks_ghost_ids[i][ghost_ids[j]] = 0
            end
        end
        for j in eachindex(put.ranks_ghost_ids[i])
            if put.ranks_ghost_ids[i][j]!=0
                push!(ut.ranks_mirror_ids[i],put.ranks_ghost_ids[i][j])
            end
        end
    end
    for i in 1:MPI.Comm_size(MPI.COMM_WORLD)
        if isempty(ut.ranks_mirror_ids[i])
            ut.mirror_datas[i] = Float64[]    
        else
            ps_datas = ut.mirror_ps_datas[ut.ranks_mirror_ids[i]]
            target_cells = ut.target_ps_datas[ut.ranks_mirror_ids[i]]
            num = sum([x.vs_data.vs_num for x in target_cells])
            ut.mirror_datas[i] = Vector{Float64}(undef,num*NDF)
            index = 1
            for j in eachindex(ps_datas)
                v = ps_datas[j].vs_data
                tv = target_cells[j].vs_data
                tdf = @views ut.mirror_datas[i][index:index+NDF*tv.vs_num-1]
                vs_project!(v.df,v.level,tdf,tv.level,v)
                index += NDF*tv.vs_num
            end
        end
        if isempty(ut.ranks_sn_ids[i])
            ut.ghost_datas[i] = Float64[]
        else
            solid_neighbors = ut.solid_neighbors[ut.ranks_sn_ids[i]]
            num = sum([x.vs_data.vs_num for x in solid_neighbors])
            ghost_id = sortperm(ut.ghost_ids[i])
            ut.ghost_datas[i] = Vector{Float64}(undef,num*NDF)
            index = 1
            for j in ghost_id
                sn = solid_neighbors[j]
                v = sn.vs_data
                sn.upwind2nd_df=reshape(@view(ut.ghost_datas[i][index:index+NDF*v.vs_num-1]),:,NDF)
                index += NDF*v.vs_num
            end
        end
    end
    amr.field.immersed_boundary.ut = ut
end
function update_mirror_data!(ut::Upwind2nd_Transport,::AMR{DIM,NDF}) where{DIM,NDF}
    for i in 1:MPI.Comm_size(MPI.COMM_WORLD)
        if !isempty(ut.ranks_mirror_ids[i])
            ps_datas = ut.mirror_ps_datas[ut.ranks_mirror_ids[i]]
            target_cells = ut.target_ps_datas[ut.ranks_mirror_ids[i]]
            index = 1
            for j in eachindex(ps_datas)
                v = ps_datas[j].vs_data
                tv = target_cells[j].vs_data
                tdf = @views ut.mirror_datas[i][index:index+NDF*tv.vs_num-1]
                vs_project!(v.df,v.level,tdf,tv.level,v)
                index += NDF*tv.vs_num
            end
        end
    end
end
function data_exchange!(ut::Upwind2nd_Transport)
    mirror_datas = ut.mirror_datas
    ghost_datas = ut.ghost_datas
    current_rank = MPI.Comm_rank(MPI.COMM_WORLD)
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(mirror_datas)
        isempty(mirror_datas[i])&&continue
        sreq = MPI.Isend(mirror_datas[i],MPI.COMM_WORLD;dest = i-1,tag = COMM_DATA_TAG+current_rank)
        push!(reqs,sreq)
    end
    for i in eachindex(ghost_datas)
        isempty(ghost_datas[i])&&continue
        rreq = MPI.Irecv!(ghost_datas[i],MPI.COMM_WORLD;source = i-1,tag = COMM_DATA_TAG+i-1)
        push!(reqs,rreq)
    end
    MPI.Waitall(reqs)
end
function data_exchange!(mirror_datas::Vector{T},ghost_datas::Vector{T}) where{T<:AbstractArray}
    current_rank = MPI.Comm_rank(MPI.COMM_WORLD)
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(mirror_datas)
        isempty(mirror_datas[i])&&continue
        sreq = MPI.Isend(mirror_datas[i],MPI.COMM_WORLD;dest = i-1,tag = COMM_DATA_TAG+current_rank)
        push!(reqs,sreq)
    end
    for i in eachindex(ghost_datas)
        isempty(ghost_datas[i])&&continue
        rreq = MPI.Irecv!(ghost_datas[i],MPI.COMM_WORLD;source = i-1,tag = COMM_DATA_TAG+i-1)
        push!(reqs,rreq)
    end
    MPI.Waitall(reqs)
end
function update_upwind2nd!(amr)
    update_local_upwind2nd!(amr.field.immersed_boundary.ut)
    update_mirror_data!(amr.field.immersed_boundary.ut,amr)
    data_exchange!(amr.field.immersed_boundary.ut)
end
function update_local_upwind2nd!(ut::Upwind2nd_Transport)
    for i in eachindex(ut.local_upwind2nd)
        upwind2nd = ut.local_upwind2nd[i]
        solid_neighbor = ut.local_solid_neighbor[i]
        df2 = upwind2nd.vs_data.df
        level2 = upwind2nd.vs_data.level
        vs_project!(df2,level2,solid_neighbor.upwind2nd_df,solid_neighbor.vs_data.level,upwind2nd.vs_data)
    end
end
function update_solid_neighbor_2nd!(::AbstractFluxType,ps_data::PS_Data{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF},amr::AMR) where{DIM,NDF}
    global_data = amr.global_data;ib = global_data.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid;opp = faceid%2==0 ? faceid-1 : faceid+1
    dir = get_dir(faceid);rot = get_rot(faceid)
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    Θ = heaviside.(vn)
    upwind1st = ps_data.neighbor.data[opp][1]
    # sL = vs_data.sdf[:,:,dir]
    upwind2nd_extrapolate!(ps_data,solid_neighbor,upwind1st,aux_df,Θ,faceid)
    # ib_df = vs_data.df+(2*aux_point[dir]-solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])*(@views vs_data.sdf[:,:,dir])
    dxL = 0.5*(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]);dxR = solid_neighbor.midpoint[dir]-aux_point[dir]
    ib_df = vs_data.df+(aux_point[dir]-dxL-ps_data.midpoint[dir])*(@views vs_data.sdf[:,:,dir])
    # ib_df = upwind1st.vs_data.df
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
    Mu_L,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
    ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
    aux_prim = IB_prim(ib,aux_point,ρw)
    F = aux_prim[1]*M
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            aux_df[i,:] .= @views F[i,:]
        end
    end
    cvc_correction!(aux_df,F,solid_neighbor,amr)
    # @. solid_neighbor.vs_data.df = 2.0*aux_df-ib_df
    @. solid_neighbor.vs_data.df =aux_df+(aux_df-ib_df)/dxL*dxR
    # @. solid_neighbor.vs_data.df = (aux_df+(aux_df-ib_df)/(aux_point[dir]-upwind1st.midpoint[dir])*
    #     (solid_neighbor.midpoint[dir]-aux_point[dir]))*Θ+(1-Θ)*(vs_data.df+(@views vs_data.sdf[:,:,dir])*
    #     (solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]))
    # @. vs_data.sdf[:,:,dir] = vanleer(sL,(solid_neighbor.vs_data.df-vs_data.df)/(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir]))
    # @. vs_data.sdf[:,:,dir] = (solid_neighbor.vs_data.df-vs_data.df)/(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])
    for i in axes(vs_data.df,1)
        if Θ[i]==1.0||vs_data.midpoint[dir]*rot>0.
            @. vs_data.sdf[i,:,dir] = @views vanleer(vs_data.sdf[i,:,dir],(solid_neighbor.vs_data.df[i,:]-vs_data.df[i,:])/dxR)
        end
    end
    solid_neighbor.w = calc_w0(vs_data.midpoint,solid_neighbor.vs_data.df,vs_data.weight,global_data)
    solid_neighbor.sw[:,dir] .= (solid_neighbor.w-ps_data.w)./(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])
end
function update_solid_neighbor_2nd!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        update_solid_neighbor_2nd!(amr.global_data.config.solver.flux,ps_data,ps_data.neighbor.data[i][1],amr)
    end
end
function update_solid_neighbor_2nd!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    update_upwind2nd!(amr)
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
            update_solid_neighbor_2nd!(ps_data,amr)
        end
    end
end
function upwind2nd_extrapolate!(ps_data::PS_Data{DIM,NDF},solid_neighbor,upwind1st,aux_df,Θ,faceid) where{DIM,NDF}
    dir = get_dir(faceid);rot = get_rot(faceid)
    vs_data = ps_data.vs_data;vn = @views vs_data.midpoint[:,dir]
    df = vs_data.df;level = vs_data.level;sdf = @views vs_data.sdf[:,:,dir]
    vs_data1 = upwind1st.vs_data;df1 = vs_data1.df;sdf1 = @views vs_data1.sdf[:,:,dir]
    level1 = vs_data1.level
    df2 = solid_neighbor.upwind2nd_df
    dx = ps_data.midpoint[dir]-upwind1st.midpoint[dir]
    flag = 0.;index = 1
    for i in axes(df,1)
        if vn[i]*rot>0.0&&Θ[i]==0.
            for j in axes(df,2)
                sdf[i,j] = (solid_neighbor.ex_df[i,j]-df[i,j])/dx
            end
            if level[i]==level1[index]
                index+=1
            elseif level[i]<level1[index]
                while flag != 1.0
                    flag += 1/2^(DIM*(level1[index]-level[i]))
                    index += 1
                end
                flag = 0.
            else
                flag += 1/2^(DIM*(level[i]-level1[index]))
                if flag == 1.
                    index += 1
                    flag = 0.
                end
            end
        else
            if level[i]==level1[index]
                for j in axes(df,2)
                    sdf[i,j] = (df[i,j]-2.0*df1[index,j]+df2[i,j])/dx+sdf1[index,j]
                    # s = (df[i,j]-2.0*df1[index,j]+df2[i,j])/dx+sdf1[index,j]
                    # r = s/(sdf[i,j]+EPS)
                    # ϕ = (r+abs(r))/(1+abs(r))
                    # sdf[i,j] = sdf[i,j]+ϕ*(s-sdf[i,j])
                    # sdf[i,j] = vanleer(s,sdf[i,j])
                end
                index += 1
            elseif level[i]<level1[index]
                sdf[i,:] .= @views (df[i,:]+df2[i,:])/dx
                while flag != 1.0
                    for j in axes(df,2)
                        sdf[i,j] += (sdf1[index,j]-2.0/dx*df1[index,j])/2^(DIM*(level1[index]-level[i]))
                    end
                    flag += 1/2^(DIM*(level1[index]-level[i]))
                    index+=1
                end
                flag = 0.
            else
                for j in axes(df,2)
                    sdf[i,j] = (df[i,j]-2.0*df1[index,j]+df2[i,j])/dx+sdf1[index,j]
                end
                flag += 1/2^(DIM*(level[i]-level1[index]))
                if flag == 1.
                    index += 1
                    flag = 0.
                end
            end
        end
    end
    @. aux_df = df + (solid_neighbor.aux_point[dir]-ps_data.midpoint[dir])*sdf
end