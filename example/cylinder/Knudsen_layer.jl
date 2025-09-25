using KitAMR,Plots,Printf
include("cylinder_udf.jl")
vs_result = KitAMR.load_vs_result("result2025-09-24_16-44")
ps_y = 0.025
ids = findall(x->x.ps_midpoint[2]≈ps_y,vs_result)
st_vs = vs_result[ids]
x = [x.ps_midpoint[1] for x in st_vs]
vs_midpoints = st_vs[1].midpoint
du = 20/60;
c = 0.5
cy = 2.5 #|3.5~6.5|
vs_id = findfirst(x->x[1]≈-c*du&&x[2]≈cy*du,eachrow(vs_midpoints))
h = [x.df[vs_id,1] for x in st_vs]
plt = plot(x,h,label="u=-"*string(@sprintf "%.2f" c*du)*",v="*string(@sprintf "%.2f" cy*du))
for i in 1:5
    ci = c+i
    vs_id = findfirst(x->x[1]≈-ci*du&&x[2]≈cy*du,eachrow(vs_midpoints))
    h = [x.df[vs_id,1] for x in st_vs]
    plot!(plt,x,h,label="u=-"*string(@sprintf "%.2f" ci*du)*",v="*string(@sprintf "%.2f" cy*du))
end
xlabel!("x")
ylabel!("h")
title!("y="*string(@sprintf "%.3f" ps_y))
display(plt)
