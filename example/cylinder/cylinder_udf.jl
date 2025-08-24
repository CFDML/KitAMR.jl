using LinearAlgebra
function shock_wave_region(midpoint,ds,global_data,level)
    if midpoint[1]>-5.0&&midpoint[1]<5.0&&midpoint[2]<5.0&&midpoint[2]>-5.0&&√sum(midpoint.^2)>1.0&&level<5
        return true
    end
    return false
end
function cylinder_buffer_IC(midpoint::Vector{Float64})
    r = norm(midpoint)
    Ma = 5.0
    Tw = 1.0
    R = 1.0
    l = 1.0 # buffer length
    if r>R+l
        return [1.0,Ma*√(5/6),0.,1.0]
    else
        return [1.0,(r-R)*Ma*√(5/6),0.,Tw-(r-R)/l*(Tw-1.0)]
    end
end