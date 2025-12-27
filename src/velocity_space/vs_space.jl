function init_vs(prim::AbstractVector,global_data::Global_Data)
    quadrature = global_data.config.quadrature
    init_vs(prim,global_data,quadrature)
end
function init_vs(vs_data::VS_Data)
    return deepcopy(vs_data)
end
function init_vs(prim::AbstractVector,global_data::Global_Data{DIM,NDF},quadrature::Vector{Float64}) where{DIM,NDF}
    trees_num = global_data.config.vs_trees_num
    vs_num = reduce(*, trees_num)
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
    # midpoints = collect(Base.Iterators.product(ranges...))
    # @inbounds @simd for i in eachindex(midpoints)
    #     midpoint[i, :] .= midpoints[i]
    # end
    midpoints = Base.Iterators.product(ranges...)
    index = 1
    for point in midpoints
        midpoint[index,:] .= point
        index += 1
    end
    df = discrete_maxwell(midpoint, prim, global_data)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    return VS_Data{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
end

function init_vs(prim::AbstractVector,global_data::Global_Data{2,NDF},quadrature::Gauss_Hermite{NP};kwargs...) where{NDF,NP}
    DIM = 2
    trees_num = global_data.config.vs_trees_num
    vs_num = reduce(*, trees_num)
    index = 1;midpoint = Matrix{Float64}(undef,vs_num,DIM);weight = Vector{Float64}(undef,vs_num)
    vcoords = quadrature.vcoords;weights = quadrature.weights
    for i in 1:NP
        for j in 1:NP
            midpoint[index,1] = vcoords[i];midpoint[index,2] = vcoords[j]
            weight[index] = weights[i]*exp(midpoint[index,1]^2)*weights[j]*exp(midpoint[index,2]^2)
            index += 1
        end
    end 
    df = haskey(kwargs,:df) ? kwargs[:df] : discrete_maxwell(midpoint, prim, global_data)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    return VS_Data{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
end