using LinearAlgebra
function shock_wave_region(midpoint,ds,global_data,level)
    refine_level = 4
    if midpoint[1]>-1.5&&midpoint[1]<1.5&&midpoint[2]<1.5&&midpoint[2]>-1.5&&√sum(midpoint.^2)>0.5&&midpoint[3]>-1.5&&midpoint[3]<1.5&&level<refine_level-2
        return true
    elseif midpoint[1]>-1.0&&midpoint[1]<1.0&&midpoint[2]<1.0&&midpoint[2]>-1.0&&√sum(midpoint.^2)>0.5&&midpoint[3]>-1.0&&midpoint[3]<1.0&&level<refine_level-1
        return true
    end
    return false
end
function sphere_buffer_IC(midpoint::Vector{Float64})
    r = norm(midpoint)
    Ma = 3.834
    Tw = 1.0+(5/3-1)*0.5*Ma^2
    if r>1.0
        return [1.0,Ma*√(5/6),0.,0.,1.0]
    else
        return [1.0,(r-0.5)/(1.0-0.5)*Ma*√(5/6),0.,0.,Tw-(r-0.5)*(Tw-1.0)]
    end
end
function vs_refine_region(midpoint;kwargs...)
    level = kwargs[:level]
    du = kwargs[:du]
    refine_level = 3
    Ma = 3.834
    Ts = 1.0+(5/3-1)*0.5*Ma^2
    U_av = [0.5*Ma*√(5/6),0.,0.]
    sn = sign(U_av-midpoint)
    midpoint_new = midpoint+0.25*du.*sn
    flag1 = level<refine_level-1;flag2 = level<refine_level
    for i in 1:3
        flag1 = flag1&&midpoint_new[i]>U_av[i]-2.0*sqrt(Ts)&&midpoint_new[i]<U_av[i]+2.0*sqrt(Ts)
        flag2 = flag2&&midpoint_new[i]>U_av[i]-sqrt(Ts)&&midpoint_new[i]<U_av[i]+sqrt(Ts)
    end
    return flag1||flag2
end