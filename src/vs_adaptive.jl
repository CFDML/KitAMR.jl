function vs_refine!(DVM_data::DVM_Data)
    trees = DVM_data.trees
    global_data = DVM_data.global_data
    vs_refine!(trees, global_data)
end
function vs_refine!(trees::Trees, global_data::Global_Data)
    ds = zeros(DIM)
    for i = 1:DIM
        ds[i] =
            (global_data.quadrature[2*i] - global_data.quadrature[2*i-1]) /
            global_data.vs_trees_num[i]
    end
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            vs_data = ps_data.vs_data
            U = ps_data.prim[2:1+DIM]
            lnmidpoint = reshape(vs_data.midpoint, :)
            lndf = reshape(vs_data.df, :)
            lnsdf = reshape(vs_data.sdf, :)
            lnflux = reshape(vs_data.flux, :)
            index = 1
            midpoint_index = Vector{Int}(undef, DIM)
            df_index = Vector{Int}(undef, NDF)
            while index < vs_data.vs_num + 1
                @inbounds for i = 1:DIM
                    midpoint_index[i] = (i - 1) * (vs_data.vs_num) + index
                end
                @inbounds for i = 1:NDF
                    df_index[i] = (i - 1) * (vs_data.vs_num) + index
                end
                midpoint = @view(lnmidpoint[midpoint_index])
                df = @view(lndf[df_index])
                if vs_data.level[index] < DVM_VS_MAXLEVEL &&
                   vs_refine_flag(ps_data.w, U, midpoint, df, vs_data.weight[index])
                    midpoint_new = midpoint_refine(midpoint, vs_data.level[index], ds)
                    df_new = df_refine(midpoint, midpoint_new, df)
                    vs_data.vs_num += 2^DIM - 1
                    level_refine_replace!(vs_data.level, index)
                    weight_refine_replace!(vs_data.weight, index)
                    midpoint_refine_replace!(
                        lnmidpoint,
                        midpoint_new,
                        vs_data.vs_num,
                        index,
                    )
                    df_refine_replace!(lndf, df_new, vs_data.vs_num, index)
                    sdf_refine_replace!(lnsdf)
                    flux_refine_replace!(lnflux)
                    index += 2^DIM - 1
                end
                index += 1
            end
            vs_data.midpoint = reshape(lnmidpoint, vs_data.vs_num, DIM)
            vs_data.df = reshape(lndf, vs_data.vs_num, NDF)
            vs_data.sdf = reshape(lnsdf, vs_data.vs_num, NDF, DIM)
            vs_data.flux = reshape(lnflux, vs_data.vs_num, NDF)
        end
    end
end
function midpoint_refine_replace!(lnmidpoint::AV, midpoint_new::AM, vs_num::Int, index::Int)
    for i = 1:DIM
        deleteat!(lnmidpoint, (i - 1) * (vs_num) + index)
        for j = 2^DIM:-1:1
            insert!(lnmidpoint, (i - 1) * (vs_num) + index, midpoint_new[j, i])
        end
    end
end
function weight_refine_replace!(weight::AV, index::Int)
    weight_new = popat!(weight, index) / 2^DIM
    for _ = 1:2^DIM
        insert!(weight, index, weight_new)
    end
end
function level_refine_replace!(level::AV, index::Int)
    level_new = popat!(level, index) + 1
    for _ = 1:2^DIM
        insert!(level, index, level_new)
    end
end
function df_refine_replace!(lndf::AV, df_new::AM, vs_num::Int, index::Int)
    for i = 1:NDF
        deleteat!(lndf, (i - 1) * (vs_num) + index)
        for j = 2^DIM:-1:1
            insert!(lndf, (i - 1) * (vs_num) + index, df_new[j, i])
        end
    end
end
function sdf_refine_replace!(lnsdf::AV)
    append!(lnsdf, zeros((2^DIM - 1) * NDF * DIM))
end
function flux_refine_replace!(lnflux::AV)
    append!(lnflux, zeros((2^DIM - 1) * NDF))
end
function vs_refine_flag(w::AV, U::AV, midpoint::AV, df::AV, weight::Float64)
    0.5 * sum((U - midpoint) .^ 2) * df[1] * weight /
    (w[4] / w[1] - 0.5 * w[1] * sum((U) .^ 2)) > 0.0001 ? true : false
end
function midpoint_refine(midpoint::AV, level::Int, ds::AV)
    midpoint_new = Matrix{Float64}(undef, 2^DIM, DIM)
    ds_new = ds / 2^(level + 1)
    for i = 1:2^DIM
        midpoint_new[i, :] .= @. midpoint + 0.5 * ds_new * rmt[i]
    end
    return midpoint_new
end
# df_refine(midpoint::AV,midpoint_new::AM,df::AV) = refine_moment_conserve(midpoint,midpoint_new)*df'
df_refine(::AV, ::AM, df::AV) = ones(2^DIM) * df'
function refine_moment_conserve(midpoint::AV, midpoint_new::AM)
    α = pinv(make_A(midpoint_new), rtol = 1e-8) * make_b(midpoint, midpoint_new)
    pushfirst!(α, 4 - sum(α))
    α
end

# function refine_moment_conserve(midpoint::AV,midpoint_new::AM)
#     model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
#     @variable(model, α[1:2^DIM-1])
#     A = make_A(midpoint_new)
#     b = make_b(midpoint,midpoint_new)
#     @objective(model, Min, sum((A * α - b).^2))
#     @constraints(model,begin
#         α .>= 0
#         sum(α)<=1.    
#     end)
#     optimize!(model)
#     α = value.(α)
#     pushfirst!(α,1-sum(α))
#     α
# end

function make_A(midpoint_new::AM)
    A = Matrix{Float64}(undef, 2^DIM - 1, 2^DIM - 1)
    A[1:DIM, :] .= transpose(@view(midpoint_new[2:end, :]))
    for i = 1:2^DIM-1
        A[end, i] = sum(@view(midpoint_new[i+1, :]) .^ 2)
    end
    for i = 1:DIM
        A[i, :] .-= midpoint_new[1, i]
    end
    A[end, :] .-= sum(@view(midpoint_new[1, :]) .^ 2)
    A
end
function make_b(midpoint::AV, midpoint_new::AM)
    b = Vector{Float64}(undef, 2^DIM - 1)
    b[1:DIM] .= midpoint .- @view(midpoint_new[1, :])
    b[end] = sum(midpoint .^ 2) .- sum(@view(midpoint_new[1, :]) .^ 2)
    b *= 2.0^DIM
    b
end

function vs_coarsen!(DVM_data::DVM_Data)
    trees = DVM_data.trees
    global_data = DVM_data.global_data
    ds = zeros(DIM)
    for i = 1:DIM
        ds[i] =
            (global_data.quadrature[2*i] - global_data.quadrature[2*i-1]) /
            global_data.vs_trees_num[i]
    end
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            vs_data = ps_data.vs_data
            U = ps_data.prim[2:1+DIM]
            lnmidpoint = reshape(vs_data.midpoint, :)
            lndf = reshape(vs_data.df, :)
            lnsdf = reshape(vs_data.sdf, :)
            lnflux = reshape(vs_data.flux, :)
            index = 1
            flag = zeros(DVM_VS_MAXLEVEL)
            midpoint_index = Matrix{Int}(undef, 2^DIM, DIM)
            df_index = Matrix{Int}(undef, 2^DIM, NDF)
            while index < vs_data.vs_num + 1
                first_level = vs_data.level[index]
                if first_level > 0
                    @inbounds for i = 1:DIM
                        midpoint_index[:, i] .=
                            (i-1)*(vs_data.vs_num)+index:(i-1)*(vs_data.vs_num)+index+2^DIM-1
                    end
                    @inbounds for i = 1:NDF
                        df_index[:, i] .=
                            (i-1)*(vs_data.vs_num)+index:(i-1)*(vs_data.vs_num)+index+2^DIM-1
                    end
                    midpoint = @view(lnmidpoint[midpoint_index])
                    df = @view(lndf[df_index])
                    if flag[first_level] % 1 == 0 &&
                       !any(x -> x > first_level, @view(vs_data.level[index+1:index+2^DIM]))
                        if vs_coarsen_flag(
                            ps_data.w,
                            U,
                            midpoint,
                            df,
                            @view(vs_data.weight[index:index+2^DIM-1])
                        )
                            midpoint_new =
                                midpoint_coarsen(@view(midpoint[1, :]), first_level, ds)
                            df_new = df_coarsen(df)
                            vs_data.vs_num -= 2^DIM - 1
                            level_coarsen_replace!(vs_data.level, index)
                            weight_coarsen_replace!(vs_data.weight, index)
                            midpoint_coarsen_replace!(
                                lnmidpoint,
                                midpoint_new,
                                vs_data.vs_num,
                                index,
                            )
                            df_coarsen_replace!(lndf, df_new, vs_data.vs_num, index)
                            sdf_coarsen_replace!(lnsdf)
                            flux_coarsen_replace!(lnflux)
                        else
                            index += 2^DIM - 1
                        end
                        if first_level > 1
                            for i = 1:first_level-1
                                flag[i] += 1 / 2^(DIM * (first_level - i))
                            end
                        end
                    else
                        for i = 1:first_level
                            flag[i] += 1 / 2^(DIM * (first_level - i + 1))
                        end
                    end
                end
                index += 1
            end
            vs_data.sdf = reshape(lnsdf, vs_data.vs_num, NDF, DIM)
            vs_data.flux = reshape(lnflux, vs_data.vs_num, NDF)
            vs_data.midpoint = reshape(lnmidpoint, vs_data.vs_num, DIM)
            vs_data.df = reshape(lndf, vs_data.vs_num, NDF)
        end
    end
end
function vs_coarsen_flag(w::AV, U::AV, midpoint::AM, df::AM, weight::AV)
    for i = 1:2^DIM
        !vs_coarsen_flag(w, U, @view(midpoint[i, :]), @view(df[i, :]), weight[i]) &&
            return false
    end
    return true
end
function vs_coarsen_flag(w::AV, U::AV, midpoint::AV, df::AV, weight::Float64)
    0.5 * sum((U - midpoint) .^ 2) * df[1] * weight /
    (w[4] / w[1] - 0.5 * w[1] * sum((U) .^ 2)) < 0.0001 / 2^DIM ? true : false
end
function midpoint_coarsen(midpoint::AV, level::Int, ds::AV)
    ds_new = ds / 2^level
    @. midpoint - 0.5 * ds_new * rmt[1]
end
function df_coarsen(df::AM)
    df_new = zeros(NDF)
    for i = 1:2^DIM
        df_new += @view(df[i, :])
    end
    df_new /= 2^DIM
    return df_new
end
function level_coarsen_replace!(level::AV, index::Int)
    for _ = 1:2^DIM-1
        deleteat!(level, index)
    end
    level[index] -= 1
end
function weight_coarsen_replace!(weight::AV, index::Int)
    for _ = 1:2^DIM-1
        deleteat!(weight, index)
    end
    weight[index] *= 2^DIM
end
function midpoint_coarsen_replace!(
    lnmidpoint::AV,
    midpoint_new::AV,
    vs_num::Int,
    index::Int,
)
    for i = 1:DIM
        for _ = 1:2^DIM-1
            deleteat!(lnmidpoint, (i - 1) * (vs_num) + index)
        end
        lnmidpoint[(i-1)*(vs_num)+index] = midpoint_new[i]
    end
end
function df_coarsen_replace!(lndf::AV, df_new::AV, vs_num::Int, index::Int)
    for i = 1:DIM
        for _ = 1:2^DIM-1
            deleteat!(lndf, (i - 1) * (vs_num) + index)
        end
        lndf[(i-1)*(vs_num)+index] = df_new[i]
    end
end
function sdf_coarsen_replace!(lnsdf::AV)
    deleteat!(lnsdf, 1:(2^DIM-1)*NDF*DIM)
end
function flux_coarsen_replace!(lnflux::AV)
    deleteat!(lnflux, 1:(2^DIM-1)*NDF)
end

function pre_vs_refine!(trees::Trees, global_data::Global_Data)
    for _ = 1:DVM_VS_MAXLEVEL
        vs_refine!(trees, global_data)
    end
end
