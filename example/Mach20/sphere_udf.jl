using LinearAlgebra
function shock_wave_region(ps_data,level,ka)
    midpoint = ps_data.midpoint;ds = ps_data.ds
    if midpoint[1]+0.5ds[1]>-2.0&&midpoint[1]-0.5ds[1]<2.0&&midpoint[2]-0.5ds[2]<2.0&&midpoint[2]+0.5ds[2]>-2.0&&midpoint[3]+0.5ds[3]>-2.0&&midpoint[3]-0.5ds[3]<2.0
        return true
    end
    return false
end
function sphere_buffer_IC(midpoint::Vector{Float64},::KInfo)
    r = norm(midpoint)
    Ma = 20.0
    Tw = 1.0+(5/3-1)*0.5*Ma^2
    if r>1.0
        return [1.0,Ma*√(5/6),0.,0.,1.0]
    else
        return [1.0,(r-0.5)/(1.0-0.5)*Ma*√(5/6),0.,0.,1.0/(Tw-(r-0.5)*(Tw-1.0))]
    end
end
function vs_comparison(ps_data,ka)
    midpoint = ps_data.midpoint
    if ps_data.vs_data.vs_num==ka.kinfo.status.max_vs_num
        return true
    end
    if (midpoint[1]>0.47-0.0625&&midpoint[1]<0.47+0.0625&&midpoint[2]>-0.59-0.0625&&midpoint[2]<-0.59+0.0625&&midpoint[3]>
            -0.0625&&midpoint[3]<eps())||(midpoint[1]>-0.5-0.0625&&midpoint[1]<-0.5&&midpoint[2]>
                eps()&&midpoint[2]<0.0625&&midpoint[3]>-0.0625&&midpoint[3]<eps())
        return true
    end
    return false
end