function ps_refine_flag(
    ps_data::PS_Data{DIM},
    amr::KitAMR_Data{DIM},
    qp::PW_pxest_quadrant_t,
) where{DIM}
    global_data = amr.global_data
    if ps_data.bound_enc!=0||domain_flag(global_data,ps_data.midpoint,ps_data.ds)
        return Cint(1)
    end
    qp.level[]>global_data.config.solver.AMR_DYNAMIC_PS_MAXLEVEL-1&&return Cint(0)
    global_data.config.user_defined.static_ps_refine_flag(ps_data.midpoint,ps_data.ds,global_data,qp.level[]) && return Cint(1)
    dflag = global_data.config.user_defined.dynamic_ps_refine_flag(;ps_data)
    !dflag&&return Cint(0)
    agrad = zeros(DIM+2)
    gradmax = global_data.status.gradmax
    for j in axes(ps_data.sw,1)
        for k in axes(ps_data.sw,2)
            agrad[j] = max(abs(ps_data.sw[j,k]),agrad[j])
        end
    end
    rgrad = maximum(agrad./gradmax)
    if rgrad > 2.0^((qp.level[] - global_data.config.solver.AMR_DYNAMIC_PS_MAXLEVEL)) * ADAPT_COEFFI_PS # 2 for second-order scheme
        flag = Cint(1)
    else
        flag = Cint(0)
    end
    flag
end

function ps_coarsen_flag(ps_datas::Vector{PS_Data}, levels::Vector{Int}, amr::KitAMR_Data{DIM,NDF}) where{DIM,NDF}
    global_data = amr.global_data
    levels[1]>global_data.config.solver.AMR_DYNAMIC_PS_MAXLEVEL&&return Cint(0)
    agrad = zeros(DIM+2)
    gradmax = global_data.status.gradmax
    for i = 1:2^DIM
        ps_data = ps_datas[i]
        (ps_data.bound_enc!=0||domain_flag(global_data,ps_data.midpoint,ps_data.ds)) && return Cint(0)
        global_data.config.user_defined.static_ps_refine_flag(ps_data.midpoint,ps_data.ds,global_data,levels[i]-1) && return Cint(0)
        for j in axes(ps_data.sw,1)
            for k in axes(ps_data.sw,2)
                agrad[j] = max(abs(ps_data.sw[j,k]),agrad[j])
            end
        end
        rgrad = maximum(agrad./gradmax)
        if rgrad > 2.0^((levels[i]-1 - global_data.config.solver.AMR_DYNAMIC_PS_MAXLEVEL)) * ADAPT_COEFFI_PS # 2 for second-order scheme
            return Cint(0)
        end
        agrad.=0.
    end
    return Cint(1)
end