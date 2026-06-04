
function initialize_vs_data(prim::AbstractVector,kinfo::KInfo)
    quadrature = kinfo.config.quadrature
    initialize_vs_data(prim,kinfo,quadrature)
end
function initialize_vs_data(vs_data::VsData)
    return deepcopy(vs_data)
end
function initialize_vs_data(prim::AbstractVector,kinfo::KInfo{DIM,NDF},quadrature::Vector{Float64}) where{DIM,NDF}
    trees_num = kinfo.config.vs_trees_num
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
    midpoints = Base.Iterators.product(ranges...)
    index = 1
    for point in midpoints
        midpoint[index,:] .= point
        index += 1
    end
    df = discrete_maxwell(midpoint, prim, kinfo)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    vs_data = VsData{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
    for _ in 1:kinfo.config.solver.AMR_VS_MAXLEVEL
        initial_vs_adaptive_mesh_refinement!(prim,vs_data,kinfo)
    end
    vs_data.df = discrete_maxwell(vs_data.midpoint, prim, kinfo)
    return vs_data
end

function initialize_vs_data(prim::AbstractVector,kinfo::KInfo{2,NDF},quadrature::Gauss_Hermite{NP};kwargs...) where{NDF,NP}
    DIM = 2
    trees_num = kinfo.config.vs_trees_num
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
    df = haskey(kwargs,:df) ? kwargs[:df] : discrete_maxwell(midpoint, prim, kinfo)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    return VsData{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
end

function initialize_vs_data(prim::AbstractVector,kinfo::KInfo{3,NDF},quadrature::Gauss_Hermite{NP};kwargs...) where{NDF,NP}
    DIM = 3
    trees_num = kinfo.config.vs_trees_num
    vs_num = reduce(*, trees_num)
    index = 1;midpoint = Matrix{Float64}(undef,vs_num,DIM);weight = Vector{Float64}(undef,vs_num)
    vcoords = quadrature.vcoords;weights = quadrature.weights
    for i in 1:NP
        for j in 1:NP
            for k in 1:NP
                midpoint[index,1] = vcoords[i];midpoint[index,2] = vcoords[j];midpoint[index,3] = vcoords[k]
                weight[index] = weights[i]*exp(midpoint[index,1]^2)*weights[j]*exp(midpoint[index,2]^2)*weights[k]*exp(midpoint[index,3]^2)
                index += 1
            end
        end
    end
    df = haskey(kwargs,:df) ? kwargs[:df] : discrete_maxwell(midpoint, prim, kinfo)
    sdf = zeros(vs_num, NDF, DIM)
    flux = zeros(vs_num, NDF)
    return VsData{DIM,NDF}(vs_num, zeros(Int, vs_num), weight, midpoint, df, sdf, flux)
end