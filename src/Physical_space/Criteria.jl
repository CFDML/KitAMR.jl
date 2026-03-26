"""
$(TYPEDSIGNATURES)
"""
function ps_refine_flag(
    ps_data::PsData{DIM},
    ka::KA{DIM},
    qp::PW_pxest_quadrant_t,
) where{DIM}
    kinfo = ka.kinfo
    if ps_data.bound_enc!=0||domain_flag(kinfo,ps_data.midpoint,ps_data.ds)
        return Cint(1)
    end
    qp.level[]>kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL-1&&return Cint(0)
    kinfo.config.user_defined.static_ps_refine_flag(ps_data.midpoint,ps_data.ds,kinfo,qp.level[]) && return Cint(1)
    dflag = kinfo.config.user_defined.dynamic_ps_refine_flag(ps_data,qp.level[],ka)
    !dflag&&return Cint(0)
    agrad = zeros(DIM+2)
    gradmax = kinfo.status.gradmax
    for j in axes(ps_data.sw,1)
        for k in axes(ps_data.sw,2)
            agrad[j] = max(abs(ps_data.sw[j,k]),agrad[j])
        end
    end
    rgrad = maximum(agrad./gradmax)
    if rgrad > 2.0^(2*(qp.level[] - kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL)) * kinfo.config.solver.ADAPT_COEFFI_PS # 2 for second-order scheme
        flag = Cint(1)
    else
        flag = Cint(0)
    end
    flag
end
"""
$(TYPEDSIGNATURES)
"""
function ps_coarsen_flag(ps_datas::Vector{PsData}, levels::Vector{Int}, ka::KA{DIM,NDF}) where{DIM,NDF}
    kinfo = ka.kinfo
    levels[1]>kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL&&return Cint(0)
    agrad = zeros(DIM+2)
    gradmax = kinfo.status.gradmax
    for i = 1:2^DIM
        ps_data = ps_datas[i]
        (ps_data.bound_enc!=0||domain_flag(kinfo,ps_data.midpoint,ps_data.ds)) && return Cint(0)
        kinfo.config.user_defined.static_ps_refine_flag(ps_data.midpoint,ps_data.ds,kinfo,levels[i]-1) && return Cint(0)
        for j in axes(ps_data.sw,1)
            for k in axes(ps_data.sw,2)
                agrad[j] = max(abs(ps_data.sw[j,k]),agrad[j])
            end
        end
        rgrad = maximum(agrad./gradmax)
        if rgrad > 2.0^(2*(levels[i]-1 - kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL)) * kinfo.config.solver.ADAPT_COEFFI_PS # 2 for second-order scheme
            return Cint(0)
        end
        agrad.=0.
    end
    return Cint(1)
end