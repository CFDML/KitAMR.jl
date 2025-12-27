using LinearAlgebra
function shock_wave_region(midpoint,ds,global_data,level)
    if midpoint[1]>-5.0&&midpoint[1]<5.0&&midpoint[2]<5.0&&midpoint[2]>-5.0&&√sum(midpoint.^2)>1.0&&level<4
        return true
    end
    return false
end
function amr_region(;ps_data,kwargs...)
    midpoint = ps_data.midpoint
    if midpoint[1]>-5.0&&midpoint[1]<5.0&&midpoint[2]<5.0&&midpoint[2]>-5.0&&√sum(midpoint.^2)>1.0
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
function vs_output_flag(;ps_data,kwargs...)
    midpoint = ps_data.midpoint;ds = ps_data.ds
    tps = Vector{Vector{Float64}}(undef,5)
    tps[1] = [-1.4589554070669333,0.05709202822871219];tps[2] = [-1.3786385660534994,0.05709202822871219]
    tps[3] = [1.7311,3.3467];tps[4] = [-0.9037743267386753,0.45412464859673063]
    tps[5] = [0.9361902761835246,0.9787846953422692]
    for i in eachindex(tps)
        if midpoint[1]-0.5ds[1]<tps[i][1]&&midpoint[1]+0.5ds[1]>tps[i][1]&&midpoint[2]-0.5ds[2]<tps[i][2]&&midpoint[2]+0.5ds[2]>tps[i][2]
            return i,true        
        end
    end
    return 0,false
end