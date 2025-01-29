"""
Adaptation flag corresponding to the velocity cell's contribution to the macro quantities. Specifically, the mass and the energy are considered.
"""
function contribution_refine_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::Global_Data{DIM,2}) where{DIM}
    max(abs(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight^2) /
    (w[end] / w[1] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight^2/w[1]) > ADAPT_COEFFI_VS ? true : false
end
function contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, global_data::Global_Data{DIM,2}) where{DIM}
    for i = 1:2^DIM
        !contribution_coarsen_flag(w, U, @view(midpoint[i, :]), @view(df[i, :]), weight[i], global_data) &&
            return false
    end
    return true
end
function contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::Global_Data{DIM,2}) where{DIM}
    max(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight^2 /
    (w[end] / w[1] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight^2/w[1]) < 1e-2* ADAPT_COEFFI_VS / 2^DIM
end
