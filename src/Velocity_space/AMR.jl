function pre_vs_refine!(trees::PS_Trees{DIM,NDF}, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    vr = vs_resolution(trees,global_data)
    !isa(global_data.config.quadrature,Vector)&&return nothing
    ds = [(global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
        global_data.config.vs_trees_num[i] for i in 1:DIM]
    vs_refine_udf = global_data.config.user_defined.static_vs_refine_flag
    for _ = 1:global_data.config.solver.AMR_VS_MAXLEVEL
        for i in eachindex(trees.data)
            for j in eachindex(trees.data[i])
                ps_data = trees.data[i][j]
                isa(ps_data,InsideSolidData) && continue
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
                        (   vs_refine_udf==null_udf ? 
                            contribution_refine_flag(ps_data.w, U, midpoint, df, vs_data.weight[index], vr, global_data) : 
                            vs_refine_udf(midpoint;level = vs_data.level[index],du = ds./2^vs_data.level[index])
                        )
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
    return nothing
end
function criterion_df(ps_data::AbstractPsData{DIM,NDF}) where{DIM,NDF}
    vs_data = ps_data.vs_data;df = vs_data.df;sdf = vs_data.sdf
    ds = ps_data.ds
    cdf = similar(vs_data.df)
    ddfi = zeros(DIM)
    for j in axes(cdf,2)
        for i in axes(cdf,1)
            for k in 1:DIM
                ddfi[k] = abs(sdf[i,j,k]*ds[k])
            end
            cdf[i,j] = df[i,j]+maximum(ddfi)
        end
    end
    return cdf
end
function vs_refine!(va_data::Velocity_Adaptive_Data, amr::KitAMR_Data{DIM,NDF}) where{DIM,NDF}
    trees = amr.field.trees;global_data = amr.global_data
    !isa(global_data.config.quadrature,Vector)&&return nothing
    ds = [(global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
    global_data.config.vs_trees_num[i] for i in 1:DIM]
    va_flags = va_data.va_flags
    id = 0
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            cdf = criterion_df(ps_data)
            vs_data = ps_data.vs_data
            U = ps_data.prim[2:1+DIM]
            lnmidpoint = reshape(vs_data.midpoint, :)
            lndf = reshape(vs_data.df, :)
            lncdf = reshape(cdf,:)
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
                df = @view(lndf[df_index]);cdf = @view(lncdf[df_index])
                if vs_data.level[index] < global_data.config.solver.AMR_VS_MAXLEVEL &&
                    contribution_refine_flag(ps_data.w, U, midpoint, cdf, vs_data.weight[index], va_data.vr,global_data)
                    midpoint_new = midpoint_refine(DIM,midpoint, vs_data.level[index], ds)
                    df_new = df_refine(DIM,midpoint, midpoint_new, df)
                    cdf_new = df_refine(DIM,midpoint,midpoint_new,df)
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
                    df_refine_replace!(DIM,NDF,lncdf,cdf_new,vs_data.vs_num,index)
                    sdf_refine_replace!(DIM,NDF,lnsdf)
                    flux_refine_replace!(DIM,NDF,lnflux)
                    index += 2^DIM - 1
                    !va_flags[id]&&(va_flags[id] = true)
                end
                index += 1
            end
            vs_data.midpoint = reshape(lnmidpoint, vs_data.vs_num, DIM)
            vs_data.df = reshape(lndf, vs_data.vs_num, NDF)
            vs_data.sdf = reshape(lnsdf, vs_data.vs_num, NDF, DIM)
            vs_data.flux = reshape(lnflux, vs_data.vs_num, NDF)
        end
    end
    return nothing
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

function vs_coarsen!(va_data::Velocity_Adaptive_Data,amr::KitAMR_Data{DIM,NDF})where{DIM,NDF}
    trees = amr.field.trees;global_data = amr.global_data
    ds = [(global_data.config.quadrature[2*i] - global_data.config.quadrature[2*i-1]) /
    global_data.config.vs_trees_num[i] for i in 1:DIM]
    va_flags = va_data.va_flags
    id = 0
    flag = zeros(global_data.config.solver.AMR_VS_MAXLEVEL)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            cdf = criterion_df(ps_data)
            vs_data = ps_data.vs_data
            U = ps_data.prim[2:1+DIM]
            lnmidpoint = reshape(vs_data.midpoint, :)
            lndf = reshape(vs_data.df, :)
            lncdf = reshape(cdf,:)
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
                        df = @view(lndf[df_index]);cdf = @view(lncdf[df_index])
                        if contribution_coarsen_flag(
                            ps_data.w,
                            U,
                            midpoint,
                            cdf,
                            @view(vs_data.weight[index:index+2^DIM-1]),
                            va_data.vr,
                            global_data
                        )
                            midpoint_new =
                                midpoint_coarsen(DIM,@view(midpoint[1, :]), first_level, ds);
                            df_new = df_coarsen(DIM,NDF,df)
                            cdf_new = df_coarsen(DIM,NDF,cdf)
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
                            df_coarsen_replace!(DIM,NDF,lncdf,cdf_new,vs_data.vs_num,index)
                            sdf_coarsen_replace!(DIM,NDF,lnsdf)
                            flux_coarsen_replace!(DIM,NDF,lnflux)
                            !va_flags[id]&&(va_flags[id] = true)
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
    return nothing
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


function vs_conserved_correction!(va_data::Velocity_Adaptive_Data,amr)
    trees = amr.field.trees
    global_data = amr.global_data
    va_flags = va_data.va_flags
    id = 0
    @inbounds for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id+=1;!va_flags[id]&&continue
            ps_data = trees.data[i][j]
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            vs_data = ps_data.vs_data
            F_c = discrete_maxwell(vs_data.midpoint, ps_data.prim, global_data)
            w = calc_w0(ps_data)
            prim = get_prim(w,global_data)
            F = discrete_maxwell(vs_data.midpoint, prim, global_data)
            vs_data.df .+= F_c-F
        end
    end
end

function vs_resolution(amr)
    trees = amr.field.trees
    vs_resolution(trees,amr.global_data)
end
function vs_resolution(trees::PS_Trees,global_data)
    density_res = 0.;energy_res = 0.
    @inbounds for tree in trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            density_res_i,energy_res_i = vs_resolution(ps_data,global_data)
            density_res = max(density_res_i,density_res)
            energy_res = max(energy_res_i,energy_res)
        end
    end
    density_res = MPI.Allreduce(density_res, (x,y)->max(x,y), MPI.COMM_WORLD)
    energy_res = MPI.Allreduce(energy_res, (x,y)->max(x,y), MPI.COMM_WORLD)
    return Velocity_Resolution(density_res,energy_res)
end
function vs_resolution(ps_data::PS_Data{DIM,1},global_data) where{DIM}
    vs_data = ps_data.vs_data
    U = ps_data.prim[2:DIM+1]
    density_max = maximum(vs_data.df)
    c2 = @views [sum((U-u).^2) for u in eachrow(vs_data.midpoint)]
    energy_max = 0.5*maximum(vs_data.df.*c2)
    vs_trees_num = global_data.config.vs_trees_num
    quadrature = global_data.config.quadrature
    du = [quadrature[2*i]-quadrature[2*i-1] for i in 1:DIM]
    weight = reduce(*,du)/reduce(*,vs_trees_num)/2^(DIM*(global_data.config.solver.AMR_VS_MAXLEVEL))
    return density_max*weight,energy_max*weight
end

function vs_adaptive_mesh_refinement!(amr)
    vr = vs_resolution(amr)
    fp = PointerWrapper(amr.global_data.forest.p4est)
    va_flags = zeros(Bool,fp.local_num_quadrants[])
    va_data = Velocity_Adaptive_Data(vr,va_flags)
    vs_refine!(va_data,amr)
    vs_coarsen!(va_data,amr)
    vs_conserved_correction!(va_data,amr)
end