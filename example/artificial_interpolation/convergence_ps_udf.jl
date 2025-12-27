using LinearAlgebra

function rotate_bc(;intersect_point,ib)
    r = intersect_point-ib.center
    ω = 1.0
    u = -ω*r[2]
    v = ω*r[1]
    return [1.,u,v,1.]
end
function artificial(midpoint)
    return [1.0+0.1*sin(π*(midpoint[1]+midpoint[2])),0.4*sin(π*(midpoint[1]+midpoint[2])),0.4*cos(π*(midpoint[1]+midpoint[2])),1.0/(1.0+0.1*cos(π*(midpoint[1]+midpoint[2])))]
end