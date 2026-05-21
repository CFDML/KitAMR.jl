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
function save_boundary_result!(ib::AbstractBoundary,ps_data,solid_neighbor::SolidNeighbor{DIM,NDF},boundary_results,ka::KA{DIM,NDF};dir_path="") where{DIM,NDF}
    kinfo = ka.kinfo;ib = kinfo.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid)
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    ib_point = aux_point+ps_data.midpoint-solid_cell.midpoint
    ib_df = image_df(ps_data,ib_point,ka)
    Θ = heaviside.(vn)
    ssdf = @views solid_neighbor.vs_data.sdf[:,:,dir]
    boundary_slope!(ssdf,vs_data.level,s_vs_data.level,ib_df,vs_data.df,
        s_vs_data.df,ib_point[dir]-ps_data.midpoint[dir],ps_data.midpoint[dir]-solid_neighbor.midpoint[dir],ka)
    dxL = aux_point[dir]-ib_point[dir]
    for j in axes(ssdf,2)
        for i in axes(ssdf,1)
            # positivity-preserving
            ssdf[i,j] = min(abs((ib_df[i,j]-eps())/(ssdf[i,j]*dxL+eps())),1.0)*ssdf[i,j]
        end
    end
    @. aux_df = ib_df+ssdf*dxL
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc;intersect_point=aux_point,ib);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,kinfo)
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
    cvc_correction!(aux_df,M,solid_neighbor,ka)
    aux_w = calc_w0(vs_data.midpoint,aux_df,vs_data.weight,kinfo)
    aux_prim = get_prim(aux_w,kinfo)
    aux_qf = calc_qf(vs_data.midpoint,aux_df,vs_data.weight,aux_prim,kinfo)
    aux_p = calc_pressure(vs_data.midpoint,aux_df,vs_data.weight,kinfo)
    push!(boundary_results[ps_data.bound_enc].midpoints,aux_point)
    push!(boundary_results[ps_data.bound_enc].normal,n)
    push!(boundary_results[ps_data.bound_enc].ps_solutions,Boundary_PS_Solution(aux_prim,aux_qf,aux_p))
    if NDF==2
        @suppress write_vs_VTK(aux_df,vs_data,ka,dir_path*"/"*string(ps_data.midpoint)*string(n),["h","b"],fieldvalues_fn)
    elseif NDF==1
        @suppress write_vs_VTK(aux_df,vs_data,ka,dir_path*"/"*string(ps_data.midpoint)*string(n),["df"],fieldvalues_fn)
    end
end
function save_boundary_result!(ib::AbstractBoundary,ps_data::PsData{DIM,NDF},boundary_results::Vector{Boundary_Solution},ka::KA{DIM,NDF};dir_path="") where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        save_boundary_result!(ib,ps_data,ps_data.neighbor.data[i][1],boundary_results,ka;dir_path)
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
function solid_cell_flag(boundaries::Vector{AbstractBoundary},midpoint::AbstractVector,ds::AbstractVector,kinfo::KInfo)
    for boundary in boundaries
        solid_cell_flag(boundary,midpoint,ds,kinfo,solid_flag(boundary,midpoint))&& return true
    end
    return false
end
function vs_extrapolate!(df::AbstractMatrix{Float64},sdf::AbstractArray{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},dx::Vector{Float64},weight::AbstractVector{Float64},::KA{DIM,NDF}) where{DIM,NDF}
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
function update_solid_cell!(::Type{T},ps_data::PsData{DIM,NDF},fluid_cells::Vector,ka::KA{DIM,NDF}) where{DIM,NDF,T<:Union{DVM,CAIDVM,UGKS}}
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
        vs_extrapolate!(fdf,fsdf,f_vs_data.level,vs_data.df,vs_data.level,dx,weight_i,ka)
    end
    ps_data.w = calc_w0(ps_data)
    ps_data.prim = get_prim(ps_data,ka.kinfo)
end
function update_solid_cell!(ps_data::PsData{DIM,NDF},ka::KA{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = [ps_data.neighbor.data[i][1] for i in fluid_dirs]
    update_solid_cell!(ka.kinfo.config.solver.flux,ps_data,fluid_cells,ka)
end



function vs_translate!(df::AbstractMatrix{Float64},level::AbstractVector{Int8},dft::AbstractMatrix{Float64},levelt::Vector{Int8},weight::AbstractVector{Float64},::KA{DIM,NDF}) where{DIM,NDF}
    j = 1;flag = 0.0
    @inbounds for i in axes(dft,1)
        if levelt[i] == level[j]
            @. dft[i, :] += @views df[j, :]*weight[i]
            j += 1
        elseif levelt[i] < level[j]
            while flag != 1.0
                @. dft[i, :] += @views df[j, :]/ 2^(DIM * (level[j] - levelt[i]))*weight[i]
                flag += 1 / 2^(DIM * (level[j] - levelt[i]))
                j += 1
            end
            flag = 0.0
        else
            @. dft[i, :] += @views df[j,:]*weight[i]
            flag += 1 / 2^(DIM * (levelt[i] - level[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
"""
$(SIGNATURES)
Filter the solid cell by fluid cells to guarantee the positivity of density and temperature.
"""
function filter_solid_cell!(ps_data::PsData{DIM,NDF},ka::KA{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = [ps_data.neighbor.data[i][1] for i in fluid_dirs]
    filter_solid_cell!(ka.kinfo.config.solver.flux,ps_data,fluid_cells,ka)
end
function filter_solid_cell!(::Type{T},ps_data::PsData{DIM,NDF},fluid_cells::Vector,ka::KA{DIM,NDF}) where{DIM,NDF,T<:Union{DVM,CAIDVM,UGKS}}
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
        fdf = f_vs_data.df
        vs_translate!(fdf,f_vs_data.level,vs_data.df,vs_data.level,weight_i,ka)
    end
    ps_data.w = calc_w0(ps_data)
    ps_data.prim = get_prim(ps_data,ka.kinfo)
end

"""
$(TYPEDSIGNATURES)
Update `solid_cell` in [`SolidNeighbor`](@ref) with interpolation.
"""
function update_solid_cell!(ka::KA)
    for ib in ka.kdata.field.immersed_boundaries
        for ps_data in ib.solid_cells
            update_solid_cell!(ps_data,ka)
        end
    end
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_solid_neighbor!(ka::KA)
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            initialize_solid_neighbor!(ps_data,ka)
        end
    end
end
function initialize_cutted_velocity_cell(n::Vector{Float64},vs_data::VsData{2},ka::KA{2,NDF}) where{NDF}
    any(x->abs(x)<1e-6,n)&&(return CuttedVelocityCells(Int[],copy(vs_data.weight),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Float64[],Float64[]))
    kinfo = ka.kinfo
    du = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i in 1:2]
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
function initialize_cutted_velocity_cell(n::AbstractVector,vs_data::VsData{3},ka::KA{3,NDF}) where{NDF}
    l = findall(x->abs(x)<1e-6,n)
    length(l)>1 &&
        (return CuttedVelocityCells(Int[],copy(vs_data.weight),Matrix{Float64}(undef,0,0),Matrix{Float64}(undef,0,0),Float64[],Float64[]))
    kinfo = ka.kinfo
    du = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i in 1:3]
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

"""
$(TYPEDSIGNATURES)
"""
function initialize_solid_neighbor!(ps_data::PsData{DIM,NDF},ka::KA{DIM,NDF}) where{DIM,NDF}
    solid_dirs = findall(x->!isnothing(x[1])&&x[1].bound_enc<0,ps_data.neighbor.data[1:2*DIM])
    vs_data = ps_data.vs_data
    for i in solid_dirs
        solid_cell = ps_data.neighbor.data[i][1]
        ps_data.bound_enc = -solid_cell.bound_enc
        ib = ka.kinfo.config.IB[-solid_cell.bound_enc]
        aux_point,normal = calc_intersect(ps_data.midpoint,solid_cell.midpoint,ps_data.ds,get_dir(i),ib)
        svsdata = VsData{DIM,NDF}(
            vs_data.vs_num,
            vs_data.level,
            vs_data.weight,
            vs_data.midpoint,
            copy(vs_data.df),
            zeros(vs_data.vs_num,NDF,DIM),
            zeros(vs_data.vs_num,NDF)
        )
        cvc = initialize_cutted_velocity_cell(normal,svsdata,ka) # heavy overhead
        ps_data.neighbor.data[i][1] = SolidNeighbor{DIM,NDF}(
            solid_cell.bound_enc,i,0,aux_point,normal,solid_cell,solid_cell.midpoint,solid_cell.ds,zeros(DIM+2),zeros(DIM+2),zeros(DIM+2,DIM),cvc,svsdata
        )
    end
    return nothing
end
function vs_interpolate!(f_df::AbstractMatrix,f_level::AbstractVector{Int8},fx,s_df,s_level,sx,b_df,bx,::KA{DIM,NDF}) where{DIM,NDF}
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

function cvc_correction!(aux_df,F::AbstractMatrix,solid_neighbor,ka)
    cvc = solid_neighbor.cvc
    for i in eachindex(cvc.indices)
        cvc.solid_dfs[i,:] .= @views F[cvc.indices[i],:]
        @views @. aux_df[cvc.indices[i],:] = (cvc.gas_weights[i]*cvc.gas_dfs[i,:]+cvc.solid_weights[i]*cvc.solid_dfs[i,:])/(cvc.gas_weights[i]+cvc.solid_weights[i])
    end
end
function image_df(ps_data,ip,ka::KA{DIM,NDF}) where{DIM,NDF}
    fluid_dirs = findall(x->!isnothing(x[1])&&!isa(x[1],AbstractInsideSolidData)&&x[1].bound_enc>=0,ps_data.neighbor.data)
    fluid_cells = Vector{AbstractPsData{DIM,NDF}}(undef,length(fluid_dirs)+1)
    for i in eachindex(fluid_dirs)
        fluid_cells[i] = ps_data.neighbor.data[fluid_dirs[i]][1]
    end
    fluid_cells[end]=ps_data
    image_df(ps_data,fluid_cells,ip,ka)
end
function image_df(ps_data::PsData{DIM,NDF},fluid_cells,ip,ka) where{DIM,NDF}
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
        vs_extrapolate!(fdf,fsdf,f_vs_data.level,image_df,vs_data.level,dx,weight_i,ka)
    end
    return image_df
end
function boundary_slope!(sdf::AbstractMatrix,level,level_n,sp_df,dc_df,sc_df,dxf::Float64,dxs::Float64,::KA{DIM,NDF}) where{DIM,NDF}
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
function update_solid_neighbor!(::Type{T},ps_data::PsData{DIM,NDF},solid_neighbor::SolidNeighbor{DIM,NDF},ka::KA) where{DIM,NDF,T<:AbstractFluxType}
    kinfo = ka.kinfo;ib = kinfo.config.IB[ps_data.bound_enc]
    vs_data = ps_data.vs_data
    aux_point = solid_neighbor.aux_point;n = solid_neighbor.normal
    faceid = solid_neighbor.faceid
    dir = get_dir(faceid)
    solid_cell = solid_neighbor.solid_cell;s_vs_data = solid_cell.vs_data
    vn = @views [dot(v,n) for v in eachrow(ps_data.vs_data.midpoint)]
    aux_df = zeros(vs_data.vs_num,NDF)
    ib_point = aux_point+ps_data.midpoint-solid_cell.midpoint
    ib_df = image_df(ps_data,ib_point,ka)
    Θ = heaviside.(vn)
    ssdf = @views solid_neighbor.vs_data.sdf[:,:,dir]
    boundary_slope!(ssdf,vs_data.level,s_vs_data.level,ib_df,vs_data.df,
        s_vs_data.df,ib_point[dir]-ps_data.midpoint[dir],ps_data.midpoint[dir]-solid_neighbor.midpoint[dir],ka)
    dxL = aux_point[dir]-ib_point[dir]
    for j in axes(ssdf,2)
        for i in axes(ssdf,1)
            # positivity-preserving
            ssdf[i,j] = min(abs((ib_df[i,j]-eps())/(ssdf[i,j]*dxL+eps())),1.0)*ssdf[i,j]
        end
    end
    @. aux_df = ib_df+ssdf*dxL
    cvc_gas_correction!(aux_df,solid_neighbor)
    aux_prim = get_bc(ib.bc;intersect_point=aux_point,ib);aux_prim[1] = 1.
    M = discrete_maxwell(vs_data.midpoint,aux_prim,kinfo)
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
    cvc_correction!(aux_df,M,solid_neighbor,ka)
    @. solid_neighbor.vs_data.df = aux_df # Positive part
    @. solid_neighbor.vs_data.flux = ssdf*(solid_neighbor.midpoint[dir]-aux_point[dir]) # Reconstruction perturbance
    # for i in axes(ssdf,1)
    #     for j in axes(ssdf,2)
    #         ssdf[i,j] = vs_data.sdf[i,j,dir] # Symmetric slope
    #     end
    # end
    # solid_neighbor.w = calc_w0(vs_data.midpoint,solid_neighbor.vs_data.df,vs_data.weight,kinfo)
    # solid_neighbor.sw[:,dir] .= (solid_neighbor.w-ps_data.w)./(solid_neighbor.midpoint[dir]-ps_data.midpoint[dir])
end
function update_solid_neighbor!(ps_data::PsData{DIM,NDF},ka::KA{DIM,NDF}) where{DIM,NDF}
    solid_neighbors = findall(x->isa(x[1],SolidNeighbor),ps_data.neighbor.data)
    for i in solid_neighbors
        update_solid_neighbor!(ka.kinfo.config.solver.flux,ps_data,ps_data.neighbor.data[i][1],ka)
    end
end
"""
$(TYPEDSIGNATURES)
Update [`VsData`](@ref) variables in [`SolidNeighbor`](@ref) with the immersed boundary method.
"""
function update_solid_neighbor!(ka::KA{DIM,NDF}) where{DIM,NDF}
    for ib in ka.kdata.field.immersed_boundaries
        for ps_data in ib.donor_cells
            update_solid_neighbor!(ps_data,ka)
        end
    end
end
function update_solid!(ka::KA{DIM,NDF}) where{DIM,NDF}
    initialize_solid_neighbor!(ka)
end

_is_ib_face(f::FullFace) = isa(f.there_data, SolidNeighbor)
_is_ib_face(f::HangingFace) = any(isa(d, SolidNeighbor) for d in f.there_data)
_is_ib_face(f::BackHangingFace) = any(isa(d, SolidNeighbor) for d in f.here_data)
_is_ib_face(::BoundaryFace) = false

"""
$(TYPEDSIGNATURES)
Collect donor cells, solid cells and IB-adjacent faces into
[`ImmersedBoundary`](@ref) objects stored in `ka.kdata.field.immersed_boundaries`.

Must be called after [`initialize_faces!`](@ref) and [`initialize_solid_neighbor!`](@ref).
"""
function initialize_immersed_boundaries!(ka::KA{DIM,NDF}) where{DIM,NDF}
    ibs = ImmersedBoundary{DIM,NDF}[]
    if isempty(ka.kinfo.config.IB)
        ka.kdata.field.immersed_boundaries = ibs
        return nothing
    end

    donors = AbstractPsData{DIM,NDF}[]
    solids = AbstractPsData{DIM,NDF}[]
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc > 0 && push!(donors, ps_data)
            ps_data.bound_enc < 0 && push!(solids, ps_data)
        end
    end

    ib_faces = AbstractFace[]
    for face in ka.kdata.field.faces
        _is_ib_face(face) && push!(ib_faces, face)
    end

    push!(ibs, ImmersedBoundary{DIM,NDF}(donors, solids, ib_faces))
    ka.kdata.field.immersed_boundaries = ibs
    return nothing
end