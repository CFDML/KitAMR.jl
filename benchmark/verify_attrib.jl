# =============================================================================
# verify_attrib.jl — (a) confirm Fix-1+2 took effect at the TYPE level, and
#                    (b) ATTRIBUTE the per-face allocations to make_face_vs /
#                        calc_flux / update_flux!, to locate the real bottleneck.
# The baseline benchmark showed Fix-1+2 barely changed allocations; this script
# shows WHY (which sub-call boxes) so we know what to fix next (RC-3 FaceVsData).
# Run: julia --project=. verify_attrib.jl     (single process; types are per-proc)
# =============================================================================
using KitAMR, MPI
using InteractiveUtils, BenchmarkTools, Printf
MPI.Init()

refine_left(midpoint, ds, kinfo, level) = (midpoint[1] < 0.0 && level < 1)
solver = Solver(; DIM=2, NDF=2, AMR_PS_MAXLEVEL=1, AMR_VS_MAXLEVEL=0,
                  PS_DYNAMIC_AMR=false, VS_DYNAMIC_AMR=false,
                  flux=CAIDVM, time_marching=CAIDVM_Marching)
gas = Gas(; K=0.0, Kn=0.075, ω=0.81, ωᵣ=0.81)
config = Configure(solver;
    geometry=[-0.5,0.5,-0.5,0.5], trees_num=[16,16],
    quadrature=[-5.0,5.0,-5.0,5.0], vs_trees_num=[16,16],
    IC=Uniform([1.,0.,0.,1.]),
    domain=[Domain(Maxwellian,1,[1.,0.,0.,1.]),
            Domain(Maxwellian,2,[1.,1.0*sqrt(5/6),0.,1.0]),
            Domain(Period,3),Domain(Period,4)],
    output=Output(solver), gas=gas, user_defined=UDF(; static_ps_refine_flag=refine_left))
p4est, ka = initialize(config)
KitAMR.slope!(p4est, ka)

KAT = typeof(ka)
println("typeof(ka) = ", KAT)

# ---- (a) confirm Fix-2: ka.kinfo chain is now concrete ----
getkinfo(ka) = ka.kinfo
getgas(ka)   = ka.kinfo.config.gas
getflux(ka)  = ka.kinfo.config.solver.flux
rt(f,tt) = (Base.return_types(f, tt)[1])
println("\n--- Fix-2 check (ka access chain) ---")
println("  ka.kinfo            :: ", rt(getkinfo,(KAT,)), isconcretetype(rt(getkinfo,(KAT,))) ? "  [concrete ✓]" : "  [abstract ✗]")
println("  ka.kinfo.config.gas :: ", rt(getgas,(KAT,)),   isconcretetype(rt(getgas,(KAT,)))   ? "  [concrete ✓]" : "  [abstract ✗]")
println("  solver.flux         :: ", rt(getflux,(KAT,)),  isconcretetype(rt(getflux,(KAT,)))  ? "  [concrete]"   : "  [abstract — read once via barrier]")

# ---- grab a real FullFace and its vs data ----
faces = ka.kdata.field.faces
ff = faces[findfirst(f->isa(f,KitAMR.FullFace), faces)]
here_vs, there_vs = KitAMR.make_face_vs(ff)
F = CAIDVM

# ---- (a) confirm Fix-1: inside the barrier, calc_flux is statically dispatched ----
# return type of calc_flux when scheme is the concrete type F:
crt = Base.return_types(KitAMR.calc_flux, (Type{CAIDVM}, typeof(here_vs), typeof(there_vs), typeof(ff), KAT))[1]
println("\n--- Fix-1 check / RC-3 gate ---")
println("  calc_flux(CAIDVM,...) return :: ", crt,
        crt==Tuple{Any,Any} ? "   <-- STILL Tuple{Any,Any}: calc_flux BODY is unstable (RC-3/RC-5)" : "")

# ---- (b) attribute per-face allocations to sub-calls (count + bytes) ----
function bench1(f, args...)
    b = @benchmark $f($(args)...) samples=2000 evals=1 seconds=8
    return (b.allocs, b.memory)
end
flux0    = KitAMR.calc_flux(F, here_vs, there_vs, ff, ka)  # warm
mfv      = bench1(KitAMR.make_face_vs, ff)
cf       = bench1(KitAMR.calc_flux, F, here_vs, there_vs, ff, ka)
flux,mic = KitAMR.calc_flux(F, here_vs, there_vs, ff, ka)
uf       = bench1(KitAMR.update_flux!, flux, mic, ff, here_vs.heavi)
total    = bench1(KitAMR.flux!, F, ff, ka)

println("\n--- per-FullFace allocation attribution (Fix-1+2 applied) ---")
@printf("  %-26s %10s %12s\n", "sub-call", "allocs", "bytes")
for (name,(a,m)) in (("make_face_vs(face)",mfv),
                     ("calc_flux(F,...)",cf),
                     ("update_flux!(...)",uf),
                     ("flux!(F,face,ka) TOTAL",total))
    @printf("  %-26s %10d %12d\n", name, a, m)
end

KitAMR.finalize!(p4est, ka)
MPI.Finalize()
