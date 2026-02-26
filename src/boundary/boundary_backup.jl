include("Circle.jl")
include("Vertices.jl")
include("Period.jl")
include("Triangles.jl")
get_bc(bc::AbstractVector;kwargs...) = copy(bc)
function get_bc(bc::Function;intersect_point,ib,kwargs...)
    bc(;intersect_point,ib)
end
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
            for ii in 1:NDF
                sdf[i, ii] = max(df_n[j,ii],0.)-df[i,ii]
            end
            j += 1
        elseif level[i] < level_n[j]
            for ii in 1:NDF
                sdf[i, ii] .= -df[i,ii]
            end
            while flag != 1.0
                for ii in 1:NDF
                    sdf[i, ii] += max(df_n[j, ii],0.)/ 2^(DIM * (level_n[j] - level[i]))
                end
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            for ii in 1:NDF
                sdf[i, ii] = max(df_n[j,ii],0.)-df[i,ii]
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
    sdf./=dx
end
function save_boundary_result!(ib::AbstractBoundary,ps_data,solid_neighbor::SolidNeighbor{DIM,NDF},boundary_results,amr::AMR{DIM,NDF};dir_path="") where{DIM,NDF}
    global_data = amr.global_data
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid)
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    # dxL = ps_data.ds[dir]
    # dxR = abs(solid_cell.midpoint[dir]-aux_point[dir])
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    # ib_point = aux_point+0.5*(ps_data.midpoint-solid_cell.midpoint)
    ib_point = aux_point+ps_data.midpoint-solid_cell.midpoint
    # ib_df =  @views vs_data.df+vs_data.sdf[:,:,dir]*(ib_point[dir]-ps_data.midpoint[dir])
    ib_df = image_df(ps_data,ib_point,amr)
    Θ = heaviside.(vn)
    # vs_interpolate!(ib_df,vs_data.level,ib_point[dir],s_vs_data.df,
    #     s_vs_data.level,solid_cell.midpoint[dir],aux_df,aux_point[dir],amr)
    ssdf = @views solid_neighbor.vs_data.sdf[:,:,dir]
    boundary_slope!(ssdf,vs_data.level,s_vs_data.level,ib_df,vs_data.df,
        s_vs_data.df,ib_point[dir]-ps_data.midpoint[dir],ps_data.midpoint[dir]-solid_neighbor.midpoint[dir])
    @. aux_df = ib_df+ssdf*(aux_point[dir]-ib_point[dir])
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc;intersect_point=aux_point,ib);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
    _,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
    ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
    # ρw = calc_IB_ρw(aux_point,ib,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
    aux_prim[1] = ρw
    M .*= aux_prim[1]
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            for j in axes(aux_df,2)
                aux_df[i,j] = M[i,j]
            end
        end
    end
    cvc_correction!(aux_df,M,solid_neighbor,amr)
    aux_w = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    aux_prim = get_prim(aux_w,global_data)
    aux_qf = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,aux_prim,global_data)
    aux_p = calc_pressure(vs_data.midpoint,aux_df,vs_data.weight,global_data)
    push!(boundary_results[ps_data.bound_enc].midpoints,aux_point)
    push!(boundary_results[ps_data.bound_enc].normal,n)
    push!(boundary_results[ps_data.bound_enc].ps_solutions,Boundary_PS_Solution(aux_prim,aux_qf,aux_p))
    if NDF==2
        @suppress write_vs_VTK(aux_df,vs_data,amr,dir_path*"/"*string(ps_data.midpoint)*string(n),["h","b"],fieldvalues_fn)
    elseif NDF==1
        @suppress write_vs_VTK(aux_df,vs_data,amr,dir_path*"/"*string(ps_data.midpoint)*string(n),["df"],fieldvalues_fn)
    end
end
function save_boundary_result!(ib::AbstractBoundary,ps_data::PS_Data{DIM,NDF},boundary_results::Vector{Boundary_Solution},amr::AMR{DIM,NDF};dir_path="") where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        save_boundary_result!(ib,ps_data,ps_data.neighbor.data[i][1],boundary_results,amr;dir_path)
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


function IB_flag(boundary::AbstractBoundary,aux_point::AbstractVector,midpoint::AbstractVector,ds::AbstractVector,level::Integer,global_data::Global_Data{DIM}) where{DIM}
    r = boundary.search_radius
    if level==global_data.config.solver.AMR_PS_MAXLEVEL
        return norm(midpoint-aux_point)<r
    else
        for i in 1:2^DIM
            norm(midpoint+0.25*ds.*RMT[DIM][i]-aux_point)<r && return true
        end
        return false
    end
end
function box_check(boundary::AbstractCircle,midpoint,ds) # out box?
    r = boundary.search_radius
    R = boundary.radius
    rx = norm(midpoint-boundary.center)
    rx>r+R+0.5*norm(ds)||rx<r-R-0.5*norm(ds)
end
function box_check(boundary::Vertices,midpoint,ds)
    box = boundary.box
    flag = false
    for i in eachindex(midpoint)
        flag = flag||(box[1][i]-ds[i]>midpoint[i]||box[2][i]+ds[i]<midpoint[i])
    end
    return boundary.solid&&flag
end
function IB_flag(aux_points::Vector{Vector{Vector{Float64}}},midpoint::AbstractVector,ds::AbstractVector,level::Int8,global_data)
    boundaries = global_data.config.IB
    for i in eachindex(boundaries)
        box_check(boundaries[i],midpoint,ds) && continue # Box check
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
function update_solid_cell!(::Union{DVM,CAIDVM,UGKS},ps_data::PS_Data{DIM,NDF},fluid_cells::Vector,amr::AMR{DIM,NDF}) where{DIM,NDF}
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
    ps_data.w = calc_w0(ps_data)
    ps_data.prim = get_prim(ps_data,amr.global_data)
end
function update_solid_cell!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = [ps_data.neighbor.data[i][1] for i in fluid_dirs]
    update_solid_cell!(amr.global_data.config.solver.flux,ps_data,fluid_cells,amr)
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
function initialize_cutted_velocity_cell(n::Vector{Float64},vs_data::VS_Data{2},amr::AMR{2,NDF}) where{NDF}
    any(x->abs(x)<1e-6,n)&&(return CuttedVelocityCells(Int[],copy(vs_data.weight),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Float64[],Float64[]))
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
function initialize_cutted_velocity_cell(n::AbstractVector,vs_data::VS_Data{3},amr::AMR{3,NDF}) where{NDF}
    l = findall(x->abs(x)<1e-6,n)
    length(l)>1 &&
        (return CuttedVelocityCells(Int[],copy(vs_data.weight),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Float64[],Float64[]))
    global_data = amr.global_data
    du = [(global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
        global_data.config.vs_trees_num[i] for i in 1:3]
    vertices = zeros(3,8)
    index = Int[];solid_weights = Float64[];gas_weights = Float64[]
    C = cut_cube_rotate(n)
    for i in 1:vs_data.vs_num
        ddu = du./2^(vs_data.level[i])
        midpoint = vs_data.midpoint[i,:]
        abs(dot(midpoint,n))>0.5*norm(ddu)&&continue
        for j in axes(vertices,2)
            vertices[:,j] .= 0.5*ANTIVT[3][j].*ddu+midpoint
        end
        flag,gas_weight,solid_weight = cut_cube(n,C,midpoint,ddu,vertices)
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
function initialize_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_dirs = findall(x->!isnothing(x[1])&&x[1].bound_enc<0,ps_data.neighbor.data[1:2*DIM])
    vs_data = ps_data.vs_data
    for i in solid_dirs
        solid_cell = ps_data.neighbor.data[i][1]
        ps_data.bound_enc = -solid_cell.bound_enc
        ib = amr.global_data.config.IB[-solid_cell.bound_enc]
        aux_point,normal = calc_intersect(ps_data.midpoint,solid_cell.midpoint,ps_data.ds,get_dir(i),ib)
        svsdata = VS_Data{DIM,NDF}(
            vs_data.vs_num,
            vs_data.level,
            vs_data.weight,
            vs_data.midpoint,
            copy(vs_data.df),
            zeros(vs_data.vs_num,NDF,DIM),
            Matrix{Float64}(undef,0,0)
        )
        cvc = initialize_cutted_velocity_cell(normal,svsdata,amr) # heavy overhead
        ps_data.neighbor.data[i][1] = SolidNeighbor{DIM,NDF}(
            solid_cell.bound_enc,i,0,aux_point,normal,solid_cell,solid_cell.midpoint,solid_cell.ds,zeros(DIM+2),zeros(DIM+2,DIM),cvc,svsdata
        )
    end
    return nothing
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
function cvc_gas_correction!(aux_df,solid_neighbor::SolidNeighbor{DIM,NDF}) where{DIM,NDF}
    cvc = solid_neighbor.cvc
    for i in eachindex(cvc.indices)
        cvc.gas_dfs[i,:] .= @views aux_df[cvc.indices[i],:]
    end
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
function image_df(ps_data,ip,amr::AMR{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = Vector{AbstractPsData{DIM,NDF}}(undef,length(fluid_dirs)+1)
    for i in eachindex(fluid_dirs)
        fluid_cells[i] = ps_data.neighbor.data[fluid_dirs[i]][1]
    end
    fluid_cells[end]=ps_data
    image_df(ps_data,fluid_cells,ip,amr)
end
function image_df(ps_data::PS_Data{DIM,NDF},fluid_cells,ip,amr) where{DIM,NDF}
    vs_data = ps_data.vs_data
    image_df = zeros(vs_data.vs_num,NDF)
    weights = Matrix{Float64}(undef,vs_data.vs_num,length(fluid_cells))
    for i in eachindex(fluid_cells)
        l = ip-fluid_cells[i].midpoint;l/=norm(l)
        weights[:,i] .= [max(0.,dot(u,l)/norm(u))^2 for u in eachrow(vs_data.midpoint)]
    end
    weight_i = Vector{Float64}(undef,vs_data.vs_num)
    weight_sum = sum(weights,dims=2)
    for i in eachindex(fluid_cells)
        f_vs_data = fluid_cells[i].vs_data
        for j in eachindex(weight_i)
            weight_i[j] = weight_sum[j]==0. ? 1.0/length(fluid_cells) : weights[j,i]/weight_sum[j]
        end
        fdf = f_vs_data.df;fsdf = f_vs_data.sdf;dx = ip-fluid_cells[i].midpoint
        vs_extrapolate!(fdf,fsdf,f_vs_data.level,image_df,vs_data.level,dx,weight_i,amr)
    end
    return image_df
end
function boundary_slope!(sdf::AbstractMatrix,level,level_n,sp_df,dc_df,sc_df,dxf::Float64,dxs::Float64)
    DIM = 3;NDF = size(sdf,2)
    sL = zeros(size(sdf,1),NDF)
    index = 1
    flag = 0.0
    df = dc_df
    dfn = sc_df
    @inbounds for i in axes(sL,1)
        if level[i] == level_n[index]
            for j = 1:NDF
                sL[i, j] = (df[i, j] - dfn[index, j]) / dxs
            end
            index += 1
        elseif level[i] < level_n[index]
            while flag != 1.0
                for j = 1:NDF
                    sL[i, j] +=
                        (df[i, j] - dfn[index, j]) / 2^(DIM * (level_n[index] - level[i])) /
                        dxs
                end
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index += 1
            end
            flag = 0.0
        else
            for j = 1:NDF
                sL[i, j] += (df[i, j] - dfn[index, j]) / dxs
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
    for j in axes(sdf,2)
        for i in axes(sdf,1)
            sdf[i,j] = minmod(sL[i,j],(sp_df[i,j]-dc_df[i,j])/dxf)
        end
    end
end
function update_solid_neighbor!(::AbstractFluxType,ps_data::PS_Data{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF},amr::AMR) where{DIM,NDF}
    global_data = amr.global_data;ib = global_data.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid)
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    ib_point = aux_point+ps_data.midpoint-solid_cell.midpoint
    ib_df = image_df(ps_data,ib_point,amr)
    Θ = heaviside.(vn)
    ssdf = @views solid_neighbor.vs_data.sdf[:,:,dir]
    boundary_slope!(ssdf,vs_data.level,s_vs_data.level,ib_df,vs_data.df,
        s_vs_data.df,ib_point[dir]-ps_data.midpoint[dir],ps_data.midpoint[dir]-solid_neighbor.midpoint[dir])
    @. aux_df = ib_df+ssdf*(aux_point[dir]-ib_point[dir])
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc;intersect_point=aux_point,ib);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,global_data)
    _,Mu_R = @views cvc_Mu(M[:,1],vn,Θ,solid_neighbor)
    ρw = cvc_density(aux_df,vn,Θ,solid_neighbor,Mu_R)
    aux_prim[1] = ρw
    M .*= aux_prim[1]
    for i in 1:vs_data.vs_num
        if Θ[i]==1.
            for j in axes(aux_df,2)
                aux_df[i,j] = M[i,j]
            end
        end
    end
    cvc_correction!(aux_df,M,solid_neighbor,amr)
    @. solid_neighbor.vs_data.df = aux_df+ssdf*(solid_neighbor.midpoint[dir]-aux_point[dir])    
    for i in axes(ssdf,1)
        for j in axes(ssdf,2)
            # ssdf[i,j] = vanleer(ssdf[i,j],(vs_data.df[i,j]-solid_neighbor.vs_data.df[i,j])/(ps_data.midpoint[dir]-solid_neighbor.midpoint[dir]))            
            ssdf[i,j] = vs_data.sdf[i,j,dir]
        end
    end
    solid_neighbor.w = calc_w0(vs_data.midpoint,solid_neighbor.vs_data.df,vs_data.weight,global_data)
    solid_neighbor.sw[:,dir] .= (solid_neighbor.w-ps_data.w)./(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])
end
function update_solid_neighbor!(ps_data::PS_Data{DIM,NDF},amr::AMR{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        update_solid_neighbor!(amr.global_data.config.solver.flux,ps_data,ps_data.neighbor.data[i][1],amr)
    end
end
function update_solid_neighbor!(amr::AMR{DIM,NDF};buffer_steps::Int=0,i::Int = typemax(Int64)) where{DIM,NDF}
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