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
    df = reshape(discrete_maxwell(midpoint, prim, global_data),:,1)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    vs_data = VS_Data{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
    return vs_data
end
function init_VS(vs_data::VS_Data)
    return deepcopy(vs_data)
end
