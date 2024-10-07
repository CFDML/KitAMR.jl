function update_solid_cell!(amr::AMR)
    boundary = amr.field.boundary
    for i in eachindex(boundary.solid_cells)
        update_solid_cell!(amr.global_data.config.IB[i],boundary.solid_cells[i],boundary.aux_points[i],boundary.IB_cells[i])
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
                @. temp[i, :] += @views df_n[j, :]/ 2^(DIM * (level[i] - level_n[j]))
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
function calc_aux_ρw()

end

function update_solid_cell!(circle::Circle,solidcells::SolidCells,::Vector{Vector{Float64}},IB_cells::IBCells)
    b = Vector{Float64}(undef,4)
    for i in eachindex(solidcells.ps_datas)
        ps_data = solidcells.ps_datas[i]
        aux_point = calc_intersect_point(circle,ps_data.midpoint)
        vs_data = first(IB_cells.IB_nodes[i]).vs_data
        n = (aux_point-circle.center)/circle.radius
        coeffi = [aux_point[1],aux_point[2],aux_point[1]*aux_point[2],1.0]
        aux_vs_temp = Vector{Matrix{Float64}}(undef,5)
        for i in eachindex(aux_vs_temp)
            aux_vs_temp[i] = copy(vs_data.df)
        end
        for j in 2:4
            vs_projection!(vs_data,IB_cells.IB_nodes[1][j].vs_data,aux_vs_temp[j])
        end
        for j in eachindex(aux_vs_temp[end])
            for k in 1:4
                b[k] = aux_vs_temp[k][j]
            end
            aux_vs_temp[end][j] = bilinear_coeffi_2D(IB_cells.IB_nodes[i][1].midpoint,IB_cells.IB_nodes[i][2].midpoint,IB_cells.IB_nodes[i][3].midpoint,IB_cells.IB_nodes[i][4].midpoint)\b*coeffi
        end   
        vn = [dot(@view(vs_data.midpoint[j,:]),-n) for j in axes(vs_data.midpoint,1)]

    end

end