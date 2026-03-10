using KitAMR, JLD2, MPI, BenchmarkTools
include("cylinder_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/overhead/configure_cylinder_dynamic.txt")
global_data = KitAMR.Global_Data(config)
p1 = load_object("./example/overhead/ps_data1.jld2");p2 = load_object("./example/overhead/ps_data2.jld2")
p1_origin = deepcopy(p1); p2_origin = deepcopy(p2)
p3 = deepcopy(p1); p4 = deepcopy(p2)
@btime begin # 1.696 ms (114198 allocations: 3.49 MiB)
    for i = 1:2
        vs_diff_slope!($p1,$p2,$global_data) # 32.611 μs (18 allocations: 112.83 KiB)
        vs_diff_flux!($p1,$p2,$global_data) # 1.582 ms (114132 allocations: 2.94 MiB)
    end
    vs_diff_iterate!($p1,$global_data) #  111.891 μs (48 allocations: 450.45 KiB)
    $p1.vs_data.df .= $p1_origin.vs_data.df
    $p1.w .= $p1_origin.w
end
@btime begin # 1.506 ms (111176 allocations: 3.09 MiB)
    for i = 1:2
        vs_same_slope!($p1,$p3,$global_data) #  14.050 μs (6 allocations: 112.42 KiB)
        vs_same_flux!($p1,$p3,$global_data) #  1.453 ms (111144 allocations: 2.76 MiB)
    end
    vs_same_iterate!($p1,$global_data) # 50.591 μs (26 allocations: 225.31 KiB)
    $p1.vs_data.df .= $p1_origin.vs_data.df
    $p1.w .= $p1_origin.w
end
@btime begin # 1.592 ms (117896 allocations: 3.27 MiB)
    for i = 1:2
        vs_same_slope!($p2,$p4,$global_data) #   14.771 μs (6 allocations: 118.92 KiB)
        vs_same_flux!($p2,$p4,$global_data) #   1.532 ms (117864 allocations: 2.92 MiB)
    end
    vs_same_iterate!($p2,$global_data) #  53.380 μs (26 allocations: 238.31 KiB)
    $p2.vs_data.df .= $p2_origin.vs_data.df
    $p2.w .= $p2_origin.w
end
function vs_diff_slope!(p1,p2,global_data)
    Ldata = Vector{KitAMR.NeighborQuad{2,2}}(undef,1);Rdata = Vector{KitAMR.NeighborQuad{2,2}}(undef,1)
    Ldata[1] = nothing;Rdata[1] = p2
    KitAMR.update_slope_inner!(Val(0),Val(1),p1,global_data,Ldata,Rdata,1)
end
function vs_same_slope!(p1,p3,global_data::KitAMR_Data.Global_Data{DIM,NDF}) where{DIM,NDF}
    dir = 1
    vs_num = p1.vs_data.vs_num
    sR = zeros(vs_num,NDF)
    sR.=(p3.vs_data.df-p1.vs_data.df)/p1.ds[dir]
    p1.vs_data.sdf[:,:,dir] .= sR
end
function vs_diff_flux!(p1,p2,global_data::KitAMR_Data.Global_Data{DIM,NDF}) where{DIM,NDF}
    rot = -1.;direction = 1;here_data = p1;there_data = p2
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    here_mid = vs_data.midpoint;there_mid = nvs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    nheavi = [x>0. for x in rot.*@views there_mid[:,direction]]
    @views here_vs = KitAMR.Face_VS_Data{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],
        vs_data.df[heavi,:],vs_data.sdf[heavi,:,:]
    )
    @views there_vs = KitAMR.Face_VS_Data{DIM,NDF}(
        nheavi,nvs_data.weight[nheavi],there_mid[nheavi,:],there_mid[nheavi,direction],
        nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
    )
    midpoint = 0.5*(p1.midpoint+p2.midpoint)
    Δt = global_data.status.Δt
    here_mid = here_vs.midpoint;there_mid = there_vs.midpoint
    here_vn = here_vs.vn;there_vn = there_vs.vn
    here_df = here_vs.df;there_df = there_vs.df
    here_sdf = here_vs.sdf;there_sdf = there_vs.sdf
    here_ps_mid = here_data.midpoint;there_ps_mid = there_data.midpoint
    @inbounds @views begin
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx[i,:],there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
        here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
        there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    end
    here_weight = here_vs.weight;there_weight = there_vs.weight
    flux = KitAMR.micro_to_macro(here_micro,here_mid,here_weight,here_data.vs_data)+KitAMR.micro_to_macro(there_micro,there_mid,there_weight,here_data.vs_data)
    area = KitAMR.face_area(here_data,direction)*rot
    KitAMR.update_macro_flux!(flux*area,here_data,there_data)
    KitAMR.update_micro_flux!(here_micro*area,there_micro*area,here_data,there_data,heavi)
end
function vs_same_flux!(p1,p3,global_data::KitAMR_Data.Global_Data{DIM,NDF}) where{DIM,NDF}
    rot = -1.;direction = 1;here_data = p1;there_data = p3
    vs_data = here_data.vs_data;nvs_data = there_data.vs_data
    here_mid = vs_data.midpoint;there_mid = nvs_data.midpoint
    heavi = [x<=0. for x in rot.*@views here_mid[:,direction]]
    nheavi = [x>0. for x in rot.*@views there_mid[:,direction]]
    @views here_vs = KitAMR.Face_VS_Data{DIM,NDF}(
        heavi,vs_data.weight[heavi],here_mid[heavi,:],here_mid[heavi,direction],
        vs_data.df[heavi,:],vs_data.sdf[heavi,:,:]
    )
    @views there_vs = KitAMR.Face_VS_Data{DIM,NDF}(
        nheavi,nvs_data.weight[nheavi],there_mid[nheavi,:],there_mid[nheavi,direction],
        nvs_data.df[nheavi,:],nvs_data.sdf[nheavi,:,:]
    )
    midpoint = 0.5*(p1.midpoint+p3.midpoint)
    Δt = global_data.status.Δt
    here_mid = here_vs.midpoint;there_mid = there_vs.midpoint
    here_vn = here_vs.vn;there_vn = there_vs.vn
    here_df = here_vs.df;there_df = there_vs.df
    here_sdf = here_vs.sdf;there_sdf = there_vs.sdf
    here_ps_mid = here_data.midpoint;there_ps_mid = there_data.midpoint
    @inbounds @views begin
        dx = [midpoint[j]-here_mid[i,j]*Δt-here_ps_mid[j] for i in axes(here_mid,1),j in axes(here_mid,2)]
        ndx = [midpoint[j]-there_mid[i,j]*Δt-there_ps_mid[j] for i in axes(there_mid,1),j in axes(there_mid,2)]
        df = [here_df[i,j]+dot(dx[i,:],here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx[i,:],there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
        here_micro = [df[i,j]*here_vn[i] for i in axes(df,1),j in axes(df,2)]
        there_micro = [ndf[i,j]*there_vn[i] for i in axes(ndf,1),j in axes(ndf,2)]
    end
    here_weight = here_vs.weight;there_weight = there_vs.weight
    flux = KitAMR.micro_to_macro(here_micro,here_mid,here_weight,here_data.vs_data)+KitAMR.micro_to_macro(there_micro,there_mid,there_weight,here_data.vs_data)
    area = KitAMR.face_area(here_data,direction)*rot
    KitAMR.update_macro_flux!(flux*area,here_data,there_data)
    index = 1;nindex = 1
    here_micro.*=area;there_micro.*=area
    @inbounds for i = 1:vs_data.vs_num
        if heavi[i]
            for j in 1:NDF
                vs_data.flux[i,j]+=here_micro[index,j]
                nvs_data.flux[i,j]-= here_micro[index,j]
            end
            index += 1
        else
            for j in 1:NDF
                vs_data.flux[i,j] += there_micro[nindex,j]
                nvs_data.flux[i,j] -= there_micro[nindex,j]
            end
            nindex+=1
        end
    end
end
function vs_diff_iterate!(ps_data,global_data)
    Δt = global_data.status.Δt;gas = global_data.config.gas
    vs_data = ps_data.vs_data
    area = reduce(*, ps_data.ds)
    ps_data.w .+= ps_data.flux .*Δt / area # Macroscopic update
    prim_c = KitAMR.get_prim(ps_data, global_data) # Conserved macroscopic variables
    f = vs_data.df
    f.+= Δt/area*vs_data.flux # Convection first
    w = KitAMR.calc_w0(vs_data.midpoint,f,vs_data.weight,global_data)
    prim = KitAMR.get_prim(w,global_data)
    τ = KitAMR.get_τ(prim_c, gas.μᵣ, gas.ω) # τ^{n+1}
    F_c = KitAMR.discrete_maxwell(vs_data.midpoint, prim_c, global_data)
    F = KitAMR.discrete_maxwell(vs_data.midpoint, prim, global_data)
    @. f += F_c-F # Conservation correction
    ps_data.qf .= qf = KitAMR.calc_qf(vs_data, prim_c) # Heatflux after convection
    F_c .+= KitAMR.shakhov_part(vs_data.midpoint, F_c, prim_c, qf, global_data) # g^{S,n+1}
    # Collision process
    f .*= τ/(τ+Δt)
    @. f += Δt/(τ+Δt)*F_c
    KitAMR.residual_check!(ps_data,prim_c,global_data)
    ps_data.prim .= prim_c
    ps_data.flux .= 0.0
    vs_data.flux .= 0.0
end
function vs_same_iterate!(ps_data,global_data)
    Δt = global_data.status.Δt;gas = global_data.config.gas
    vs_data = ps_data.vs_data
    area = reduce(*, ps_data.ds)
    ps_data.w .+= ps_data.flux .*Δt / area # Macroscopic update
    f = vs_data.df
    f.+= Δt/area*vs_data.flux # Convection first
    prim = KitAMR.get_prim(ps_data.w,global_data)
    τ = KitAMR.get_τ(prim, gas.μᵣ, gas.ω) # τ^{n+1}
    F = KitAMR.discrete_maxwell(vs_data.midpoint, prim, global_data)
    ps_data.qf .= qf = KitAMR.calc_qf(vs_data, prim) # Heatflux after convection
    f .*= τ/(τ+Δt)
    @. f += Δt/(τ+Δt)*F
    KitAMR.residual_check!(ps_data,prim,global_data)
    ps_data.flux .= 0.0
    vs_data.flux .= 0.0
end