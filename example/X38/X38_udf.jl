using LinearAlgebra
function amr_region(;ps_data,kwargs...)
    midpoint = ps_data.midpoint
    if midpoint[1]>-2.0&&midpoint[1]<3.0&&midpoint[2]<2.0&&midpoint[2]>-2.0&&midpoint[3]<2.0&&midpoint[3]>-2.0
        return true
    end
    return false
end
function X38_buffer_IC(midpoint::Vector{Float64};kwargs...)
    global_data = kwargs[:global_data]
    ib = global_data.config.IB[1]
    Ma = 8.
    Tw = 300/56
    _,distance = KitAMR.nn(ib.tkdt.kdt,midpoint)
    if distance > 1.0
        return [1.0,Ma*√(5/6),0.,0.,1.0]
    else
        return [1.0,distance*Ma*√(5/6),0.,0.,1.0/(1.0+(1.0-distance)*(Tw-1.0))]
    end
end
function vs_refine_region(midpoint;kwargs...)
    level = kwargs[:level]
    du = kwargs[:du]
    refine_level = 2;DIM = 3
    Ma = 3.834
    Ts = 1.0+(5/3-1)*0.5*Ma^2
    U_av = [0.5*Ma*√(5/6),0.,0.]
    sn = sign.(U_av-midpoint)
    midpoint_new = midpoint+0.25*du.*sn
    flag1 = level<refine_level-1;flag2 = level<refine_level
    for i in 1:DIM
        flag1 = flag1&&midpoint_new[i]>U_av[i]-2.0*sqrt(Ts)&&midpoint_new[i]<U_av[i]+2.0*sqrt(Ts)
        flag2 = flag2&&midpoint_new[i]>U_av[i]-1.0*sqrt(Ts)&&midpoint_new[i]<U_av[i]+1.0*sqrt(Ts)
    end
    return flag1||flag2
end