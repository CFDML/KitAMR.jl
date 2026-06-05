"""
$(TYPEDSIGNATURES)
"""
function contribution_refine_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    vr::Velocity_Resolution,
    global_data::Global_Data,
)
    return local_contribution_refine_flag(w, U, midpoint, df, weight, global_data)||global_contribution_refine_flag(
        U,
        midpoint,
        df,
        weight,
        vr,
        global_data,
    )
end
"""
$(TYPEDSIGNATURES)
"""
function contribution_coarsen_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractMatrix,
    df::AbstractMatrix,
    weight::AbstractVector,
    vr::Velocity_Resolution,
    global_data::Global_Data,
)
    return local_contribution_coarsen_flag(w, U, midpoint, df, weight, global_data)&&global_contribution_coarsen_flag(
        U,
        midpoint,
        df,
        weight,
        vr,
        global_data,
    )
end

function local_contribution_refine_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    ::Global_Data{DIM,2},
) where {DIM}
    max(
        abs(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2]) * weight) /
        (w[end] - 0.5 * w[1] * sum((U) .^ 2)),
        df[1]*weight/w[1],
    ) > 1e-2 ? true : false
end
function local_contribution_refine_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    ::Global_Data{DIM,1},
) where {DIM}
    max(
        abs(0.5 * sum((U - midpoint) .^ 2) * df[1] * weight) /
        (w[end] - 0.5 * w[1] * sum((U) .^ 2)),
        df[1]*weight/w[1],
    ) > 1e-2 ? true : false
end
function global_contribution_refine_flag(
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    vr::Velocity_Resolution,
    ::Global_Data{DIM,2},
) where {DIM}
    df[1]*weight > ADAPT_COEFFI_VS*vr.density ||
        0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2]) * weight>vr.energy*ADAPT_COEFFI_VS
end
function global_contribution_refine_flag(
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    vr::Velocity_Resolution,
    ::Global_Data{DIM,1},
) where {DIM}
    df[1]*weight > vr.density*ADAPT_COEFFI_VS ||
        0.5 * (sum((U - midpoint) .^ 2) * df[1]) * weight>vr.energy*ADAPT_COEFFI_VS
end
function local_contribution_coarsen_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractMatrix,
    df::AbstractMatrix,
    weight::AbstractVector,
    global_data::Global_Data{DIM,2},
) where {DIM}
    for i = 1:(2^DIM)
        !local_contribution_coarsen_flag(
            w,
            U,
            @view(midpoint[i, :]),
            @view(df[i, :]),
            weight[i],
            global_data,
        ) && return false
    end
    return true
end
function local_contribution_coarsen_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractMatrix,
    df::AbstractMatrix,
    weight::AbstractVector,
    global_data::Global_Data{DIM,1},
) where {DIM}
    for i = 1:(2^DIM)
        !local_contribution_coarsen_flag(
            w,
            U,
            @view(midpoint[i, :]),
            @view(df[i, :]),
            weight[i],
            global_data,
        ) && return false
    end
    return true
end
function local_contribution_coarsen_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    ::Global_Data{DIM,2},
) where {DIM}
    max(
        0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2]) * weight /
        (w[end] - 0.5 * w[1] * sum((U) .^ 2)),
        df[1]*weight/w[1],
    ) < 1e-2 / 2^DIM
end
function local_contribution_coarsen_flag(
    w::AbstractVector,
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    ::Global_Data{DIM,1},
) where {DIM}
    max(
        0.5 * (sum((U - midpoint) .^ 2) * df[1]) * weight /
        (w[end] - 0.5 * w[1] * sum((U) .^ 2)),
        df[1]*weight/w[1],
    ) < 1e-2 / 2^DIM
end

function global_contribution_coarsen_flag(
    U::AbstractVector,
    midpoint::AbstractMatrix,
    df::AbstractMatrix,
    weight::AbstractVector,
    vr::Velocity_Resolution,
    global_data::Global_Data{DIM,2},
) where {DIM}
    for i = 1:(2^DIM)
        !global_contribution_coarsen_flag(
            U,
            @view(midpoint[i, :]),
            @view(df[i, :]),
            weight[i],
            vr,
            global_data,
        ) && return false
    end
    return true
end
function global_contribution_coarsen_flag(
    U::AbstractVector,
    midpoint::AbstractMatrix,
    df::AbstractMatrix,
    weight::AbstractVector,
    vr::Velocity_Resolution,
    global_data::Global_Data{DIM,1},
) where {DIM}
    for i = 1:(2^DIM)
        !global_contribution_coarsen_flag(
            U,
            @view(midpoint[i, :]),
            @view(df[i, :]),
            weight[i],
            vr,
            global_data,
        ) && return false
    end
    return true
end
function global_contribution_coarsen_flag(
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    vr::Velocity_Resolution,
    ::Global_Data{DIM,2},
) where {DIM}
    df[1]*weight < ADAPT_COEFFI_VS*vr.density/2^DIM &&
        0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2]) * weight<vr.energy*ADAPT_COEFFI_VS/2^DIM
end
function global_contribution_coarsen_flag(
    U::AbstractVector,
    midpoint::AbstractVector,
    df::AbstractVector,
    weight::Float64,
    vr::Velocity_Resolution,
    ::Global_Data{DIM,1},
) where {DIM}
    df[1]*weight < vr.density*ADAPT_COEFFI_VS/2^DIM &&
        0.5 * (sum((U - midpoint) .^ 2) * df[1]) * weight<vr.energy*ADAPT_COEFFI_VS/2^DIM
end

function macro_estimate_refine_flag(prim::AbstractVector, U, midpoint, ds, level)
    return false
    if norm(midpoint-U)-norm(ds)/2^(level+1)<√(1/prim[end])
        return true
    else
        return false
    end
end
function macro_estimate_IB_refine_flag(prim, bc, U, midpoint, ds, level)
    Ub = bc[2:(length(midpoint)+1)]
    if (
        norm(midpoint-U)-norm(ds)/2^(level+1)<√(1/prim[end])||norm(midpoint-Ub)-norm(ds)/2^(
            level+1
        )<√(1/bc[end])
    )
        return true
    else
        return false
    end
end
