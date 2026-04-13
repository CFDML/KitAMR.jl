using LinearAlgebra
function amr_region(ps_data,level,ka) # dynamic
    midpoint = ps_data.midpoint
    ds = ps_data.ds
    L = 1.0
    if midpoint[1]+ds[1]>-2.0*L&&midpoint[1]-ds[1]<2.0*L&&midpoint[2]-ds[2]<2.0*L&&midpoint[2]+ds[2]>-2.0*L&&midpoint[3]-ds[3]<2.0*L&&midpoint[3]+ds[3]>-2.0*L
        return true
    end
    return false
end
function X38_buffer_IC(midpoint::Vector{Float64},kinfo::KInfo)
    ib = kinfo.config.IB[1]
    Ma = 8.
    T0 = 1.0
    Tw = 300/56
    _,distance = KitAMR.nn(ib.tkdt.kdt,midpoint)
    R = 0.3
    if distance > R
        return [1.0,Ma*√(5/6),0.,0.,1/T0]
    else
        return [1.0,distance*Ma*√(5/6)/R,0.,0.,1.0/(T0+(R-distance)*(Tw-T0)/R)]
    end
end