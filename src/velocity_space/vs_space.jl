function init_VS(prim::AbstractVector,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    trees_num = global_data.config.vs_trees_num
    vs_num = reduce(*, trees_num)
    quadrature = global_data.config.quadrature
    midpoint = zeros(vs_num, DIM)
    ds = zeros(DIM)
    ranges = Array{StepRangeLen}(undef, DIM)
    for i = 1:DIM
        half_len = (quadrature[2*i] - quadrature[2*i-1]) / (2 * trees_num[i])
        ds[i] = 2 * half_len
        ranges[i] =
            range(quadrature[2*i-1] + half_len, quadrature[2*i] - half_len, trees_num[i])
    end
    weight = reduce(*, ds) .* ones(vs_num)
    midpoints = collect(Base.Iterators.product(ranges...))
    @inbounds @simd for i in eachindex(midpoints)
        midpoint[i, :] .= midpoints[i]
    end
    df = discrete_maxwell(midpoint, prim, global_data)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    vs_data = VS_Data{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
    return vs_data
end
function init_VS(vs_data::VS_Data)
    return deepcopy(vs_data)
end
function vs_project(df::AbstractMatrix,level::AbstractVector,tlevel::AbstractVector,vs_data::AbstractVsData{DIM,NDF}) where{DIM,NDF}
    tdf = Matrix{Float64}(undef,vs_data.vs_num,NDF)
    vs_project!(df,level,tdf,tlevel,vs_data)
    return tdf
end
function vs_project!(df::AbstractMatrix,level::AbstractVector,tdf::AbstractVector,tlevel::AbstractVector,vs_data::AbstractVsData{DIM,NDF}) where{DIM,NDF}
    tdf = reshape(tdf,:,NDF)
    vs_project!(df,level,tdf,tlevel,vs_data)
end
function vs_project!(df::AbstractMatrix,level::AbstractVector,tdf::AbstractMatrix,tlevel::AbstractVector,::AbstractVsData{DIM,NDF}) where{DIM,NDF}
    index = 1;flag = 0.
    @inbounds for i in axes(tdf,1)
        if tlevel[i]==level[index]
            for j in 1:NDF
                tdf[i,j] = df[index,j]
            end
            index += 1
        elseif tlevel[i] < level[index]
            tdf[i,:] .= 0.
            while flag != 1.0
                for j in 1:NDF
                    tdf[i,j] += df[index,j]/2^(DIM*(level[index]-tlevel[i]))
                end
                flag += 1/2^(DIM*(level[index]-tlevel[i]))
                index += 1
            end
            flag = 0.
        else
            for j in 1:NDF
                tdf[i,j] = df[index,j]
            end
            flag += 1/2^(DIM*(tlevel[i]-level[index]))
            if flag == 1.
                index += 1
                flag = 0.
            end
        end
    end
end