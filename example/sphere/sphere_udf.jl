using LinearAlgebra
function shock_wave_region(midpoint,ds,global_data,level)
    refine_level = 5
    if midpoint[1]>-5.0&&midpoint[1]<5.0&&midpoint[2]<5.0&&midpoint[2]>-5.0&&√sum(midpoint.^2)>1.0&&midpoint[3]>-5.0&&midpoint[3]<5.0&&level<5
        return true
    end
    return false
end
function sphere_buffer_IC(midpoint::Vector{Float64})
    r = norm(midpoint)
    Ma = 2.0
    Tw = 1.0
    if r>2.0
        return [1.0,Ma*√(5/6),0.,0.,1.0]
    else
        return [1.0,(r-1.0)*Ma*√(5/6),0.,0.,Tw-(r-1.0)*(Tw-1.0)]
    end
end