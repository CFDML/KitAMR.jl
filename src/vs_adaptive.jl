function vs_refine!(amr::AMR)
    trees = amr.field.trees
    global_data = amr.global_data
    vs_refine!(trees, global_data)
end
function vs_refine!(trees::PS_Trees{DIM,NDF}, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    ds = zeros(DIM)
    for i = 1:DIM
        ds[i] =
            (global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
            global_data.config.vs_trees_num[i]
    end
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            # ps_data.bound_enc<0 && continue # solid_cell
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
                if vs_data.level[index] < global_data.config.solver.AMR_VS_MAXLEVEL &&
                   vs_refine_flag(ps_data.w, U, midpoint, df, vs_data.weight[index], global_data)
                    midpoint_new = midpoint_refine(DIM,midpoint, vs_data.level[index], ds)
                    df_new = df_refine(DIM,midpoint, midpoint_new, df)
                    vs_data.vs_num += 2^DIM - 1
                    level_refine_replace!(DIM,vs_data.level, index)
                    weight_refine_replace!(DIM,vs_data.weight, index)
                    midpoint_refine_replace!(
                        DIM,
                        lnmidpoint,
                        midpoint_new,
                        vs_data.vs_num,
                        index,
                    )
                    df_refine_replace!(DIM,NDF,lndf, df_new, vs_data.vs_num, index)
                    sdf_refine_replace!(DIM,NDF,lnsdf)
                    flux_refine_replace!(DIM,NDF,lnflux)
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
function midpoint_refine_replace!(DIM::Integer,lnmidpoint::AbstractVector, midpoint_new::AbstractMatrix, vs_num::Int, index::Int)
    for i = 1:DIM
        deleteat!(lnmidpoint, (i - 1) * (vs_num) + index)
        for j = 2^DIM:-1:1
            insert!(lnmidpoint, (i - 1) * (vs_num) + index, midpoint_new[j, i])
        end
    end
end
function weight_refine_replace!(DIM::Integer,weight::AbstractVector, index::Int)
    weight_new = popat!(weight, index) / 2^DIM
    for _ = 1:2^DIM
        insert!(weight, index, weight_new)
    end
end
function level_refine_replace!(DIM::Integer,level::AbstractVector, index::Int)
    level_new = popat!(level, index) + 1
    for _ = 1:2^DIM
        insert!(level, index, level_new)
    end
end
function df_refine_replace!(DIM::Integer,NDF::Integer,lndf::AbstractVector, df_new::AbstractMatrix, vs_num::Int, index::Int)
    for i = 1:NDF
        deleteat!(lndf, (i - 1) * (vs_num) + index)
        for j = 2^DIM:-1:1
            insert!(lndf, (i - 1) * (vs_num) + index, df_new[j, i])
        end
    end
end
function sdf_refine_replace!(DIM::Integer,NDF::Integer,lnsdf::AbstractVector)
    append!(lnsdf, zeros((2^DIM - 1) * NDF * DIM))
end
function flux_refine_replace!(DIM::Integer,NDF::Integer,lnflux::AbstractVector)
    append!(lnflux, zeros((2^DIM - 1) * NDF))
end
function vs_refine_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::Global_Data{DIM,2}) where{DIM}
    max(abs(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight^2) /
    (w[end] / w[1] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight^2/w[1]) > ADAPT_COEFFI_VS ? true : false
end
function midpoint_refine(DIM::Integer,midpoint::AbstractVector, level::Int8, ds::AbstractVector)
    midpoint_new = Matrix{Float64}(undef, 2^DIM, DIM)
    ds_new = ds / 2^(level + 1)
    for i = 1:2^DIM
        midpoint_new[i, :] .= @. midpoint + 0.5 * ds_new * RMT[DIM][i]
    end
    return midpoint_new
end
# df_refine(DIM::Integer,midpoint::AbstractVector,midpoint_new::AbstractMatrix,df::AbstractVector) = refine_moment_conserve(midpoint,midpoint_new)*df'
df_refine(DIM::Integer,::AbstractVector, ::AbstractMatrix, df::AbstractVector) = ones(2^DIM) * df'
function refine_moment_conserve(DIM::Integer,midpoint::AbstractVector, midpoint_new::AbstractMatrix)
    α = pinv(make_A(DIM,midpoint_new), rtol = 1e-8) * make_b(DIM,midpoint, midpoint_new)
    pushfirst!(α, 4 - sum(α))
    α
end

# function refine_moment_conserve(midpoint::AbstractVector,midpoint_new::AbstractMatrix)
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

function make_A(DIM::Integer,midpoint_new::AbstractMatrix)
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
function make_b(DIM::Integer,midpoint::AbstractVector, midpoint_new::AbstractMatrix)
    b = Vector{Float64}(undef, 2^DIM - 1)
    b[1:DIM] .= midpoint .- @view(midpoint_new[1, :])
    b[end] = sum(midpoint .^ 2) .- sum(@view(midpoint_new[1, :]) .^ 2)
    b *= 2.0^DIM
    b
end

function vs_coarsen!(amr::AMR{DIM,NDF})where{DIM,NDF}
    trees = amr.field.trees
    global_data = amr.global_data
    ds = zeros(DIM)
    for i = 1:DIM
        ds[i] =
            (global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
            global_data.config.vs_trees_num[i]
    end
    flag = zeros(global_data.config.solver.AMR_VS_MAXLEVEL)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            # ps_data.bound_enc<0 && continue # solid_cell
            vs_data = ps_data.vs_data
            U = ps_data.prim[2:1+DIM]
            lnmidpoint = reshape(vs_data.midpoint, :)
            lndf = reshape(vs_data.df, :)
            lnsdf = reshape(vs_data.sdf, :)
            lnflux = reshape(vs_data.flux, :)
            index = 1;flag.=0.
            midpoint_index = Matrix{Int}(undef, 2^DIM, DIM)
            df_index = Matrix{Int}(undef, 2^DIM, NDF)
            while index < vs_data.vs_num + 1
                first_level = vs_data.level[index]
                if first_level > 0
                    if flag[first_level]%1==0. &&
                       all(x -> x == first_level, @view(vs_data.level[index+1:index+2^DIM-1]))
                        for i in axes(midpoint_index,2)
                            midpoint_index[:, i] .=
                                (i-1)*(vs_data.vs_num)+index:(i-1)*(vs_data.vs_num)+index+2^DIM-1
                        end
                        for i in axes(df_index,2)
                            df_index[:, i] .=
                                (i-1)*(vs_data.vs_num)+index:(i-1)*(vs_data.vs_num)+index+2^DIM-1
                        end
                        midpoint = @view(lnmidpoint[midpoint_index])
                        df = @view(lndf[df_index])
                        if vs_coarsen_flag(
                            ps_data.w,
                            U,
                            midpoint,
                            df,
                            @view(vs_data.weight[index:index+2^DIM-1]),
                            global_data
                        )
                            midpoint_new =
                                midpoint_coarsen(DIM,@view(midpoint[1, :]), first_level, ds)
                            df_new = df_coarsen(DIM,NDF,df)
                            vs_data.vs_num -= 2^DIM - 1
                            level_coarsen_replace!(DIM,vs_data.level, index)
                            weight_coarsen_replace!(DIM,vs_data.weight, index)
                            midpoint_coarsen_replace!(
                                DIM,
                                lnmidpoint,
                                midpoint_new,
                                vs_data.vs_num,
                                index,
                            )
                            df_coarsen_replace!(DIM,NDF,lndf, df_new, vs_data.vs_num, index)
                            sdf_coarsen_replace!(DIM,NDF,lnsdf)
                            flux_coarsen_replace!(DIM,NDF,lnflux)
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
                            flag[i] += 1 / 2^(DIM * (first_level - i+1))
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
function vs_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, global_data::Global_Data{DIM,2}) where{DIM}
    for i = 1:2^DIM
        !vs_coarsen_flag(w, U, @view(midpoint[i, :]), @view(df[i, :]), weight[i], global_data) &&
            return false
    end
    return true
end
function vs_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::Global_Data{DIM,2}) where{DIM}
    max(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight^2 /
    (w[end] / w[1] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight^2/w[1]) < 1e-2* ADAPT_COEFFI_VS / 2^DIM
end
function midpoint_coarsen(DIM::Integer,midpoint::AbstractVector, level::Int8, ds::AbstractVector)
    ds_new = ds / 2^level
    @. midpoint - 0.5 * ds_new * RMT[DIM][1]
end
function df_coarsen(DIM::Integer,NDF::Integer,df::AbstractMatrix)
    df_new = zeros(NDF)
    for i = 1:2^DIM
        df_new += @view(df[i, :])
    end
    df_new /= 2^DIM
    return df_new
end
function level_coarsen_replace!(DIM::Integer,level::AbstractVector, index::Int)
    for _ = 1:2^DIM-1
        deleteat!(level, index)
    end
    level[index] -= 1
end
function weight_coarsen_replace!(DIM::Integer,weight::AbstractVector, index::Int)
    for _ = 1:2^DIM-1
        deleteat!(weight, index)
    end
    weight[index] *= 2^DIM
end
function midpoint_coarsen_replace!(
    DIM::Integer,
    lnmidpoint::AbstractVector,
    midpoint_new::AbstractVector,
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
function df_coarsen_replace!(DIM::Integer,NDF::Integer,lndf::AbstractVector, df_new::AbstractVector, vs_num::Int, index::Int)
    for i = 1:NDF
        for _ = 1:2^DIM-1
            deleteat!(lndf, (i - 1) * (vs_num) + index)
        end
        lndf[(i-1)*(vs_num)+index] = df_new[i]
    end
end
function sdf_coarsen_replace!(DIM::Integer,NDF::Integer,lnsdf::AbstractVector)
    deleteat!(lnsdf, 1:(2^DIM-1)*NDF*DIM)
end
function flux_coarsen_replace!(DIM::Integer,NDF::Integer,lnflux::AbstractVector)
    deleteat!(lnflux, 1:(2^DIM-1)*NDF)
end

function pre_vs_refine!(trees::PS_Trees{DIM,NDF}, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    for _ = 1:global_data.config.solver.AMR_VS_MAXLEVEL
        vs_refine!(trees, global_data)
    end
end
# function update_solid_vs!(amr::AMR)
#     boundary = amr.field.boundary
#     solid_cells = boundary.solid_cells
#     IB_cells = boundary.IB_cells
#     @inbounds for i in eachindex(solid_cells)
#         for j in eachindex(solid_cells[i].ps_datas)
#             ps_data = solid_cells[i].ps_datas[j]
#             vs_data = ps_data.vs_data
#             IB_vs = first(IB_cells[i].IB_nodes[j]).vs_data
#             vs_data.vs_num = IB_vs.vs_num
#         end
#     end
# end