function init_VS(prim::AV, global_data::Global_Data)
    trees_num = global_data.vs_trees_num
    vs_num = reduce(*, trees_num)
    quadrature = global_data.quadrature
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
    df = discrete_maxwell(midpoint, prim, global_data.gas.K)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    vs_data = VS_Data(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
    # neighbor_states = Vector{Vector{Vector{Int}}}(undef,2*DIM)
    # @inbounds @simd for i = 1:2*DIM
    #     neighbor_states[i] = [ones(Int,vs_num)]
    # end
    # vs_structure = VS_Structure(neighbor_states,ones(Int,vs_num),ones(Int,vs_num))
    return vs_data
end
function init_VS(vs_data::VS_Data)
    return deepcopy(vs_data)
end
