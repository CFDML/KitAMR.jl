# In gaussian weighted adaptive method, ps_data.sw and vs_data.flux[:,1] are respectively used to store prims(the aux_pointed (by fluids) one and the interpolated one) and Θ
function gaussian_weight(λ::Float64,midpoint_in::AbstractVector,midpoint_out::AbstractVector) # here, midpoint represents the thermal velocity
    exp(-λ*(sum(midpoint_in.^2)-sum(midpoint_out.^2)))
end
function df_refine(df::AbstractVector,midpoints::AbstractMatrix,midpoint::AbstractVector,U::AbstractVector,prim::AbstractVector,::PS_Data{DIM,NDF}) where{DIM,NDF}
    λ = prim[end]
    weights = [gaussian_weight(λ,x-U,midpoint-U) for x in eachrow(midpoints)]
    df_new = [weights[i]*df[j] for i in 1:2^DIM, j in 1:NDF]
    return df_new
end
function df_coarsen(df::AbstractMatrix,midpoint::AbstractVector,midpoints::AbstractMatrix,U::AbstractVector,prim::AbstractVector)
    λ = prim[end]
    weight = exp(-λ*sum((midpoint-U).^2))/sum([exp(-λ*sum((x-U).^2)) for x in eachrow(midpoints)])
    df_new = vec(sum(df;dims=1)*weight)
    return df_new
end
# function flux_coarsen_replace!(index,lnflux,::PS_Data{DIM,NDF}) where{DIM,NDF}
#     for _ = 1:2^DIM-1
#         deleteat!(lnflux, index)
#     end
#     # lnflux[index] = df_new[i]
#     for _ in 1:(NDF-1)*(2^DIM-1)
#         pop!(lnflux)
#     end
# end
# function flux_refine_replace(index,lnflux,::PS_Data{DIM,NDF}) where{DIM,NDF}
#     Θ = lnflux[index]
#     for _ = 1:2^DIM-1
#         insert!(lnflux,index,Θ)
#     end
#     for _ in 1:(NDF-1)*(2^DIM-1)
#         push!(0,lnflux)
#     end
# end