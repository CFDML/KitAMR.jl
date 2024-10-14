function RBF_coeffi(points::AbstractVector{Vector{T}}) where{T}
    N = length(points)
    dis = 0
    for i in 1:N-1
        for j in i+1:N
            dis += norm(points[i]-points[j])
        end
    end
    sigma = dis/(N*(N-1)/2)
    
    return RBF_coeffi(sigma,points),sigma
end
function RBF_coeffi(sigma::Real,points::AbstractVector{Vector{T}}) where{T}
    N = length(points)
    Φ = Matrix{T}(undef,N,N)
    for i in 1:N
        for j in 1:N
            if j==i
                Φ[i,j] = 1.0
            else
                Φ[i,j] = exp(-sigma^2*norm(points[i]-points[j])^2)
            end
        end
    end
    return Φ
end

function update_solid_cell!(amr::AMR)
    global_data = amr.global_data
    boundary = amr.field.boundary
    for i in eachindex(boundary.solid_cells)
        update_solid_cell!(global_data.config.IB[i],boundary.solid_cells[i],boundary.aux_points[i],boundary.IB_cells[i],global_data)
    end
end

function vs_projection!(vs_data::VS_Data,vs_data_n::VS_Data,temp::AbstractMatrix)
    j = 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    df_n = vs_data_n.df
    for i = 1:vs_data.vs_num
        if level[i] == level_n[j]
            @. temp[i, :] .= @views df_n[j, :]
            j += 1
        elseif level[i] < level_n[j]
            while flag != 1.0
                @. temp[i, :] += @views df_n[j, :]/ 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @. temp[i, :] = @views df_n[j,:]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end

function bilinear_coeffi_2D(p1::T,p2::T,p3::T,p4::T) where{T<:AbstractVector{Float64}}
    return [
        p1[1] p1[2] p1[1]*p1[2] 1.0;
        p2[1] p2[2] p2[1]*p2[2] 1.0;
        p3[1] p3[2] p3[1]*p3[2] 1.0;
        p4[1] p4[2] p4[1]*p4[2] 1.0
    ]
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Circle,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_2D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end

function calc_IB_ρw_2D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_2D2F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, vn, Θ)
end
function IB_prim(circle::Circle,aux_point::AbstractVector,ρw::Real)
    IB_prim(circle.bc,aux_point,ρw)
end
function IB_prim(bc::AbstractVector,::AbstractVector,ρw::Real)
    prim = copy(bc)
    prim[1] = ρw
    return prim
end
# function normalize_point(ref_point::AbstractVector,scale::Real,points::AbstractVector{T}...) where{T}
#     norm_points = Vector{Vector{T}}(undef,length(points))
#     for i in eachindex(norm_points)
#         norm_points[i] = (points[i]-ref_point)/scale
#     end
#     return norm_points
# end
function make_coeffi!(coeffi::AbstractVector,sigma::Real,interp_point::AbstractVector{T},Vpoints::AbstractVector{Vector{T}}) where{T}
    for i in eachindex(coeffi)
        coeffi[i] = exp(-sigma^2*norm(interp_point-Vpoints[i])^2)
    end
end
function update_solid_cell!(circle::Circle,solidcells::SolidCells{DIM,NDF},::Vector{Vector{Float64}},IB_cells::IBCells,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    b = Vector{Float64}(undef,4)
    for i in eachindex(solidcells.ps_datas)
        ps_data = solidcells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        image_point = 2*aux_point-ps_data.midpoint
        vs_data = first(IB_cells.IB_nodes[i]).vs_data
        n = (aux_point-circle.center)/circle.radius
        aux_vs_temp = Vector{Matrix{Float64}}(undef,6)
        dxL = norm(aux_point-image_point)
        dxR = norm(ps_data.midpoint-aux_point)
        for i in eachindex(aux_vs_temp)
            aux_vs_temp[i] = Matrix{Float64}(undef,vs_data.vs_num,NDF)
        end
        aux_vs_temp[1] .= vs_data.df
        ip_df = aux_vs_temp[end-1]
        aux_df = aux_vs_temp[end]
        for j in 2:4
            vs_projection!(vs_data,IB_cells.IB_nodes[1][j].vs_data,aux_vs_temp[j])
        end
        vn = [dot(@view(vs_data.midpoint[j,:]),n) for j in axes(vs_data.midpoint,1)]
        Θ = heaviside.(vn)
        # points = normalize_point(aux_point,dxL,IB_cells.IB_nodes[i][1].midpoint,IB_cells.IB_nodes[i][2].midpoint,IB_cells.IB_nodes[i][3].midpoint,IB_cells.IB_nodes[i][4].midpoint)
        points = [IB_cells.IB_nodes[i][j].midpoint for j in 1:4]
        A,sigma = RBF_coeffi(points)
        ip_coeffi = Vector{Float64}(undef,4);aux_coeffi = Vector{Float64}(undef,4)
        make_coeffi!(ip_coeffi,sigma,image_point,points);make_coeffi!(aux_coeffi,sigma,aux_point,points)
        # calc distribution function at intersect point
        for j in axes(aux_df,1)
            for l in 1:NDF
                for k in 1:4
                    b[k] = aux_vs_temp[k][j,l]
                end
                RBF_weights = A\b
                if Θ[j]<0
                    aux_df[j,l] = dot(RBF_weights,aux_coeffi)
                end
                ip_df[j,l] = dot(RBF_weights,ip_coeffi)
            end
        end   
        ρw = calc_IB_ρw(aux_point,circle,vs_data.midpoint,vs_data.weight,vs_data.df,vn,Θ)
        aux_prim = IB_prim(circle,aux_point,ρw)
        for j in axes(aux_df,1)
            if Θ[j]<=0
                aux_df[j,:] .= discrete_maxwell(@view(vs_data.midpoint[j,:]),aux_prim,global_data)
            end
        end
        ip_prim = Vector{Float64}(undef,4)
        for j in 1:4
            for k in 1:4
                b[k] = IB_cells.IB_nodes[i][k].prim[j]
            end
            ip_prim[j] = dot(A\b,ip_coeffi)
        end
        ps_data.prim = @. aux_prim + (aux_prim-ip_prim)*dxR/dxL
        ps_data.w = get_conserved(ps_data,global_data)
        s_vs_data = ps_data.vs_data
        s_vs_data.vs_num = vs_data.vs_num
        s_vs_data.midpoint = vs_data.midpoint
        s_vs_data.level = vs_data.level
        s_vs_data.weight = vs_data.weight
        s_vs_data.df = @. aux_df[end]+(aux_df-ip_df)*dxR/dxL
        !(size(s_vs_data.sdf)==size(vs_data.sdf)) && (s_vs_data.sdf = Array{Float64}(undef,size(vs_data.sdf)))
    end
end
function project_solid_cell_slope!(vs_data::AbstractVsData{DIM,NDF},vs_data_n::VS_Data{DIM,NDF},DIR::Integer) where{DIM,NDF}
    j = 1
    flag = 0.0
    level = vs_data.level
    sdf = @view(vs_data.sdf[:,:,DIR])
    level_n = vs_data_n.level
    sdf_n = @view(vs_data_n.sdf[:,:,DIR])
    for i = 1:vs_data.vs_num
        if level[i] == level_n[j]
            @. sdf[i, :] .= @views sdf_n[j, :]
            j += 1
        elseif level[i] < level_n[j]
            while flag != 1.0
                @. sdf[i, :] += @views sdf_n[j, :]/ 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @. sdf[i, :] = @views sdf_n[j,:]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end