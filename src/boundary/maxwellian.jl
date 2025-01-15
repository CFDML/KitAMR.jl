function WLS_coeffi(interp_point::AbstractVector{T},aux_point::AbstractVector{T},points::AbstractVector{Vector{T}},R::T) where{T} # interp_point: the point to be interpolated (image_point or aux_point); R: search radius
    # w_m is calculated according to the distance to the aux_point
    N = length(points);DIM = length(interp_point)
    w = Vector{T}(undef,N+1)
    # _,index = findmin(norm(interp_point.-point) for point in points)
    interp_point = interp_point-aux_point
    points = [point-aux_point for point in points]
    # points[1],points[index] = points[index],points[1]
    for i in 1:N
        w[i] = 0.5*(1+cos(π*norm(points[i])/R))
    end
    w[end] = 1.0
    W = Diagonal(w)
    M = WLST[DIM-1]
    N<M && @error `WLS lacks points!`
    V = Matrix{T}(undef,N+1,M)
    if DIM==2
        for i in 1:N
            V[i,1] = 1.0
            V[i,2] = points[i][1]
            V[i,3] = points[i][2]
            V[i,4] = points[i][1]^2
            V[i,5] = points[i][2]^2
            V[i,6] = points[i][1]*points[i][2]
        end
        V[end,1] = 1.0
        V[end,2] = interp_point[1]
        V[end,3] = interp_point[2]
        V[end,4] = interp_point[1]^2
        V[end,5] = interp_point[2]^2
        V[end,6] = interp_point[1]*interp_point[2]
    elseif DIM==3
        for i in 1:N
            V[i,1] = 1.0
            V[i,2] = points[i][1]
            V[i,3] = points[i][2]
            V[i,4] = points[i][3]
            V[i,5] = points[i][1]^2
            V[i,6] = points[i][2]^2
            V[i,7] = points[i][3]^2
            V[i,8] = points[i][1]*points[i][2]
            V[i,9] = points[i][1]*points[i][3]
            V[i,10] = points[i][2]*points[i][3]
        end
        V[end,1] = 1.0
        V[end,2] = interp_point[1]
        V[end,3] = interp_point[2]
        V[end,4] = interp_point[3]
        V[end,5] = interp_point[1]^2
        V[end,6] = interp_point[2]^2
        V[end,7] = interp_point[3]^2
        V[end,8] = interp_point[1]*interp_point[2]
        V[end,9] = interp_point[1]*interp_point[3]
        V[end,10] = interp_point[2]*interp_point[3]
    end
    return pinv(W*V;atol = 1e-5)*W
end
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
                Φ[i,j] = exp(-norm(points[i]-points[j])^2/sigma^2)
            end
        end
    end
    return Φ
end

function update_solid_cell!(amr::AMR)
    global_data = amr.global_data
    boundary = amr.field.boundary
    for i in eachindex(boundary.solid_cells)
        @inbounds update_solid_cell!(global_data.config.IB[i],boundary.solid_cells[i],boundary.aux_points[i],boundary.IB_cells[i],global_data)
    end
end

function vs_projection!(vs_data::AbstractVsData{DIM,NDF},vs_data_n::AbstractVsData,temp::AbstractMatrix) where{DIM,NDF}
    j = 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    df_n = vs_data_n.df
    @inbounds for i = 1:vs_data.vs_num
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
    return @inbounds [
        p1[1] p1[2] p1[1]*p1[2] 1.0;
        p2[1] p2[2] p2[1]*p2[2] 1.0;
        p3[1] p3[2] p3[1]*p3[2] 1.0;
        p4[1] p4[2] p4[1]*p4[2] 1.0
    ]
end
function bilinear_coeffi_3D(p1::T,p2::T,p3::T,p4::T,p5::T,p6::T,p7::T) where{T<:AbstractVector{Float64}}
    return @inbounds [
        p1[1] p1[2] p1[3] p1[1]*p1[2] p1[1]*p1[3] p1[2]*p1[3] 1.0;
        p2[1] p2[2] p2[3] p2[1]*p2[2] p2[1]*p2[3] p2[2]*p2[3] 1.0;
        p3[1] p3[2] p3[3] p3[1]*p3[2] p3[1]*p3[3] p3[2]*p3[3] 1.0;
        p4[1] p4[2] p4[3] p4[1]*p4[2] p4[1]*p4[3] p4[2]*p4[3] 1.0;
        p5[1] p5[2] p5[3] p5[1]*p5[2] p5[1]*p5[3] p5[2]*p5[3] 1.0;
        p6[1] p6[2] p6[3] p6[1]*p6[2] p6[1]*p6[3] p6[2]*p6[3] 1.0;
        p7[1] p7[2] p7[3] p7[1]*p7[2] p7[1]*p7[3] p7[2]*p7[3] 1.0
    ]
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Circle,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_2D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Sphere,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_3D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end

function calc_IB_ρw_2D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_2D2F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function calc_IB_ρw_3D(::AbstractVector,bc::AbstractVector,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    maxwellian_density_3D1F(@view(midpoint[:,1]),@view(midpoint[:,2]),@view(df[:,1]), bc, weight, Θ, vn)
end
function IB_prim(circle::AbstractCircle,aux_point::AbstractVector,ρw::Real)
    IB_prim(circle.bc,aux_point,ρw)
end
function IB_prim(bc::AbstractVector,::AbstractVector,ρw::Real)
    prim = copy(bc)
    prim[1] = ρw
    return prim
end
# function normalize_point(ref_point::AbstractVector,scale::Real,points::Vector{Vector{T}}) where{T}
#     norm_points = Vector{Vector{T}}(undef,length(points))
#     for i in eachindex(norm_points)
#         norm_points[i] = (points[i]-ref_point)/scale
#     end
#     return norm_points
# end
# function make_coeffi!(coeffi::AbstractVector,sigma::Real,interp_point::AbstractVector{T},Vpoints::AbstractVector{Vector{T}}) where{T}
#     for i in eachindex(coeffi)
#         coeffi[i] = exp(-norm(interp_point-Vpoints[i])^2/sigma^2)
#     end
# end
function make_coeffi!(coeffi::AbstractVector,sigma::Real,interp_point::AbstractVector{T},Vpoints::AbstractVector{Vector{T}}) where{T}
    for i in eachindex(coeffi)
        coeffi[i] = exp(-norm(interp_point-Vpoints[i])^2/sigma^2)
    end
end
function make_bilinear_coeffi_2D(image_point::AbstractVector)
    return [image_point[1],image_point[2],image_point[1]*image_point[2],1.]
end
# function update_solid_cell!(circle::Circle,solidcells::SolidCells{DIM,NDF},::Vector{Vector{Float64}},IB_cells::IBCells,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
#     # b = Vector{Float64}(undef,4)
#     for i in eachindex(solidcells.ps_datas)
#         ps_data = solidcells.ps_datas[i]
#         aux_point = calc_intersect_point(circle,ps_data.midpoint)
#         image_point = 2*aux_point-ps_data.midpoint
#         vs_data = first(IB_cells.IB_nodes[i]).vs_data
#         n = (aux_point-circle.center)/circle.radius # outer normal direction
#         M = length(IB_cells.IB_nodes[i])
#         aux_vs_temp = Vector{Matrix{Float64}}(undef,M+2)
#         dxL = norm(aux_point-image_point)
#         dxR = norm(ps_data.midpoint-aux_point)
#         for i in eachindex(aux_vs_temp)
#             aux_vs_temp[i] = Matrix{Float64}(undef,vs_data.vs_num,NDF)
#         end
#         aux_vs_temp[1] .= vs_data.df
#         ip_df = aux_vs_temp[end-1]
#         aux_df = aux_vs_temp[end]
#         for j in 2:length(IB_cells.IB_nodes[i])
#             vs_projection!(vs_data,IB_cells.IB_nodes[i][j].vs_data,aux_vs_temp[j])
#         end
#         vn = [dot(@view(vs_data.midpoint[j,:]),n) for j in axes(vs_data.midpoint,1)]
#         Θ = heaviside.(vn)
#         # points = normalize_point(aux_point,dxL,IB_cells.IB_nodes[i][1].midpoint,IB_cells.IB_nodes[i][2].midpoint,IB_cells.IB_nodes[i][3].midpoint,IB_cells.IB_nodes[i][4].midpoint)
#         points = [IB_cells.IB_nodes[i][j].midpoint for j in eachindex(IB_cells.IB_nodes[i])]
#         # A,sigma = RBF_coeffi(points)
#         # image_point interpolate by WLS
#         A_ip = WLS_coeffi(image_point,aux_point,points,circle.search_radius)
#         # ip_coeffi = Vector{Float64}(undef,4);aux_coeffi = Vector{Float64}(undef,4)
#         # make_coeffi!(ip_coeffi,sigma,image_point,points);make_coeffi!(aux_coeffi,sigma,aux_point,points)
#         # calc distribution function at intersect point
#         # _,index_ip = findmin(norm(image_point.-point) for point in points)
#         # iter = filter(i->i!=index_ip,1:M)
#         for j in axes(ip_df,1)
#             for l in 1:NDF
#                 ip_df[j,l] = (aux_vs_temp[1][j,l]-sum(A_ip[end,m]*aux_vs_temp[m][j,l] for m in 1:M))/A_ip[end,end]
#             end
#         end   
#         ip_prim = Vector{Float64}(undef,4)
#         for j in 1:4
#             ip_prim[j] = (IB_cells.IB_nodes[i][1].prim[j]-sum(A_ip[end,m]*IB_cells.IB_nodes[i][m].prim[j] for m in 1:M))/A_ip[end,end]
#         end
#         # aux_point interpolate by WLS
#         A_ap = WLS_coeffi(aux_point,aux_point,points,circle.search_radius)
#         for j in axes(aux_df,1)
#             if Θ[j]==0.
#                 for l in 1:NDF
#                     aux_df[j,l] = (aux_vs_temp[1][j,l]-sum(A_ap[end,m]*aux_vs_temp[m][j,l] for m in 1:M))/A_ap[end,end]
#                 end
#             end
#         end   
#         ρw = calc_IB_ρw(aux_point,circle,vs_data.midpoint,vs_data.weight,aux_df,vn,Θ)
#         # ρ1 = sum(vs_data.weight.*@view(vs_data.df[:,1]))
        
#         # ρw_test = calc_IB_ρw(aux_point,circle,vs_data.midpoint,vs_data.weight,vs_data.df,vn,Θ)
#         aux_prim = IB_prim(circle,aux_point,ρw)
#         for j in axes(aux_df,1)
#             if Θ[j]==1.
#                 aux_df[j,:] .= discrete_maxwell(@view(vs_data.midpoint[j,:]),aux_prim,global_data)
#             end
#         end
#         # sum1 = sum(@view(vs_data.df[:,1]))
#         # sumap = sum(@view(aux_df[:,1]))
#         # sumM = [sum(@view(aux_vs_temp[m][:,1])) for m in 2:M]
#         # @show minimum(@view(ip_df[:,1])) minimum(@view(aux_df[:,1]))
#         # dfs = [aux_vs_temp[m][340,1] for m in 1:M]
#         # if i==2
#         #     @show findmin(@views aux_df[:,1]) vs_data.midpoint[340,:] Θ[340] dfs points aux_point aux_prim ρw_test ρ1 n
#         # end
#         if ρw>1.5||ρw<0||isnan(ρw)
#             @show ρw points aux_point
#         end
#         # @show ρw aux_prim A_ip[end,end]
#         # @show aux_prim sumap sumM sum1 
#         ps_data.prim .= @. aux_prim + (aux_prim-ip_prim)*dxR/dxL
#         ps_data.w = get_conserved(ps_data,global_data)
#         s_vs_data = ps_data.vs_data
#         s_vs_data.vs_num = vs_data.vs_num
#         s_vs_data.midpoint = vs_data.midpoint
#         s_vs_data.level = vs_data.level
#         s_vs_data.weight = vs_data.weight
#         s_vs_data.df = @. aux_df[end]+(aux_df-ip_df)*dxR/dxL
#         !(size(s_vs_data.sdf)==size(vs_data.sdf)) && (s_vs_data.sdf = Array{Float64}(undef,size(vs_data.sdf)))
#     end
# end
function update_solid_cell!(circle::Circle,solidcells::SolidCells{DIM,NDF},::Vector{Vector{Float64}},IB_cells::IBCells,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    for i in eachindex(solidcells.ps_datas)
        ps_data = solidcells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        image_point = 2*aux_point-ps_data.midpoint
        # IB_vs = first(IB_cells.IB_nodes[i]).vs_data
        s_vs_data = ps_data.vs_data
        n = (aux_point-circle.center)/circle.radius # outer normal direction
        aux_vs_temp = Vector{Matrix{Float64}}(undef,6)
        dxL = norm(aux_point-image_point)
        dxR = norm(ps_data.midpoint-aux_point)
        for i in eachindex(aux_vs_temp)
            aux_vs_temp[i] = zeros(Float64,s_vs_data.vs_num,NDF)
        end
        # aux_vs_temp[1] .= IB_vs.df
        ip_df = aux_vs_temp[end-1]
        aux_df = aux_vs_temp[end]
        for j in 1:4
            vs_projection!(s_vs_data,IB_cells.IB_nodes[i][j].vs_data,aux_vs_temp[j])
        end
        vn = [dot(@view(s_vs_data.midpoint[j,:]),n) for j in axes(s_vs_data.midpoint,1)]
        Θ = heaviside.(vn)
        points = [IB_cells.IB_nodes[i][j].midpoint for j in 1:4]
        # try Ainv = inv(bilinear_coeffi_2D(points...))
        # catch
        #     @show length(points) points
        # end
        Ainv = inv(bilinear_coeffi_2D(points...))
        b = Vector{Float64}(undef,4);ip_coeffi = make_bilinear_coeffi_2D(image_point)
        @inbounds for j in axes(ip_df,1)
            for l in 1:NDF
                for k in 1:4
                    b[k] = aux_vs_temp[k][j,l]
                end
                ip_df[j,l] = dot(Ainv*b,ip_coeffi)
            end
        end   
        # aux_point interpolate by bilinear
        ap_coeffi = make_bilinear_coeffi_2D(aux_point)
        @inbounds for j in axes(aux_df,1)
            if Θ[j]==0.
                for l in 1:NDF
                    for k = 1:4
                        b[k] = aux_vs_temp[k][j,l]
                    end
                    aux_df[j,l] = dot(Ainv*b,ap_coeffi)
                end
            end
        end   
        ρw = calc_IB_ρw(aux_point,circle,s_vs_data.midpoint,s_vs_data.weight,aux_df,vn,Θ)
        aux_prim = IB_prim(circle,aux_point,ρw)
        # sdf = (aux_df-ip_df)/dxL
        for j in axes(aux_df,1)
            if Θ[j]==1.
                aux_df[j,:] .= discrete_maxwell(@view(s_vs_data.midpoint[j,:]),aux_prim,global_data)
            end
        end
        # s_vs_data = ps_data.vs_data
        # s_vs_data.vs_num = vs_data.vs_num
        # s_vs_data.midpoint = vs_data.midpoint
        # s_vs_data.level = vs_data.level
        # s_vs_data.weight = vs_data.weight
        # s_vs_data.df = @. aux_df+sdf*dxR
        s_vs_data.df = @. aux_df+(aux_df-ip_df)/dxL*dxR
        # s_vs_data.sdf
		ps_data.w = calc_w0(ps_data)
		ps_data.prim = get_prim(ps_data,global_data)
        # size(s_vs_data.sdf,1)!=vs_data.vs_num && (s_vs_data.sdf = Array{Float64}(undef,vs_data.vs_num,NDF,DIM))
    end
end
function update_solid_cell!(circle::Sphere,solidcells::SolidCells{DIM,NDF},::Vector{Vector{Float64}},IB_cells::IBCells,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    for i in eachindex(solidcells.ps_datas)
        ps_data = solidcells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        image_point = 2*aux_point-ps_data.midpoint
        s_vs_data = ps_data.vs_data
        n = (aux_point-circle.center)/circle.radius # outer normal direction
        aux_vs_temp = Vector{Matrix{Float64}}(undef,9)
        dxL = norm(aux_point-image_point)
        dxR = norm(ps_data.midpoint-aux_point)
        for i in eachindex(aux_vs_temp)
            aux_vs_temp[i] = zeros(Float64,s_vs_data.vs_num,NDF)
        end
        ip_df = aux_vs_temp[end-1]
        aux_df = aux_vs_temp[end]
        for j in 1:7
            vs_projection!(s_vs_data,IB_cells.IB_nodes[i][j].vs_data,aux_vs_temp[j])
        end
        vn = [dot(@view(s_vs_data.midpoint[j,:]),n) for j in axes(s_vs_data.midpoint,1)]
        Θ = heaviside.(vn)
        points = [IB_cells.IB_nodes[i][j].midpoint for j in 1:7]
        Ainv = inv(bilinear_coeffi_3D(points...))
        b = Vector{Float64}(undef,7);ip_coeffi = make_bilinear_coeffi_3D(image_point)
        @inbounds for j in axes(ip_df,1)
            for l in 1:NDF
                for k in 1:7
                    b[k] = aux_vs_temp[k][j,l]
                end
                ip_df[j,l] = dot(Ainv*b,ip_coeffi)
            end
        end   
        # aux_point interpolate by bilinear
        ap_coeffi = make_bilinear_coeffi_2D(aux_point)
        @inbounds for j in axes(aux_df,1)
            if Θ[j]==0.
                for l in 1:NDF
                    for k = 1:7
                        b[k] = aux_vs_temp[k][j,l]
                    end
                    aux_df[j,l] = dot(Ainv*b,ap_coeffi)
                end
            end
        end   
        ρw = calc_IB_ρw(aux_point,circle,s_vs_data.midpoint,s_vs_data.weight,aux_df,vn,Θ)
        aux_prim = IB_prim(circle,aux_point,ρw)
        for j in axes(aux_df,1)
            if Θ[j]==1.
                aux_df[j,:] .= discrete_maxwell(@view(s_vs_data.midpoint[j,:]),aux_prim,global_data)
            end
        end
        s_vs_data.df = @. aux_df+(aux_df-ip_df)/dxL*dxR
		ps_data.w = calc_w0(ps_data)
		ps_data.prim = get_prim(ps_data,global_data)
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
            @. sdf[i, :] = @views sdf_n[j, :]
            j += 1
        elseif level[i] < level_n[j]
            sdf[i, :] .= 0.
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
function calc_solid_cell_slope!(svdata::AbstractVsData{DIM,NDF},fvdata::VS_Data{DIM,NDF},smid::Vector{Float64},fmid::Vector{Float64},direction::Integer) where{DIM,NDF}
    j = 1
    flag = 0.0
    level = svdata.level
    sdf = @view(svdata.sdf[:,:,direction])
    df = svdata.df
    level_n = fvdata.level
    df_n = fvdata.df
    dx = fmid[direction]-smid[direction]
    for i in 1:svdata.vs_num
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
