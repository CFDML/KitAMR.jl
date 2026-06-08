# Type-stability diagnosis of flux!(face, ka) — Cthulhu-style manual descent via code_warntype.
# Run: mpirun -np 1 julia --project=. diag_typestab.jl
using KitAMR, MPI
using InteractiveUtils
MPI.Init()

const RANK = MPI.Comm_rank(MPI.COMM_WORLD)

# ---- refine the left half at init so the mesh has hanging faces ----
refine_left(midpoint, ds, kinfo, level) = (midpoint[1] < 0.0 && level < 1)

solver = Solver(;
    DIM = 2, NDF = 2,
    AMR_PS_MAXLEVEL = 1,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
)
gas = Gas(; K = 0.0, Kn = 0.075, ω = 0.81, ωᵣ = 0.81)
output = Output(solver)
udf = UDF(; static_ps_refine_flag = refine_left)
config = Configure(solver;
    geometry = [-0.5,0.5,-0.5,0.5],
    trees_num = [16,16],
    quadrature = [-5.0,5.0,-5.0,5.0],
    vs_trees_num = [16,16],
    IC = Uniform([1.,0.,0.,1.]),
    domain = [
        Domain(Maxwellian,1,[1.,0.,0.,1.]),
        Domain(Maxwellian,2,[1.,1.0*sqrt(5/6),0.,1.0]),
        Domain(Period,3), Domain(Period,4),
    ],
    output = output, gas = gas, user_defined = udf,
)

p4est, ka = initialize(config)
KitAMR.slope!(p4est, ka)   # populate slopes so calc_flux paths are realistic

faces = ka.kdata.field.faces
println("="^70)
println("rank $RANK : typeof(ka) = ", typeof(ka))
println("rank $RANK : #faces = ", length(faces))
# tally concrete face types actually present
let d = Dict{Any,Int}()
    for f in faces; d[typeof(f)] = get(d, typeof(f), 0) + 1; end
    for (k,v) in d; println("   present: ", k, "  x", v); end
end
println("="^70)

# ---------- helper: capture code_warntype text (no color) ----------
function warntype_str(f, tt)
    buf = IOBuffer()
    io = IOContext(buf, :color => false)
    try
        code_warntype(io, f, tt)
    catch e
        println(io, "ERROR running code_warntype: ", e)
    end
    return String(take!(buf))
end

# ---------- helper: programmatic instability scan via code_typed ----------
# Returns (return_type, list of (ssa_index, abstract_type)) for non-concrete SSA values.
function scan_unstable(f, tt)
    out = Tuple{Int,Any}[]
    local rt
    try
        ci_rt = code_typed(f, tt; optimize = true)
        isempty(ci_rt) && return (Any, out, "no method")
        ci, rt = ci_rt[1]
        sts = ci.ssavaluetypes
        if sts isa Vector
            for (i, t) in enumerate(sts)
                t === Any && (push!(out, (i, Any)); continue)
                if t isa Type && !isconcretetype(t) && t !== Union{} &&
                   !(t <: Type) && t !== Nothing
                    push!(out, (i, t))
                end
            end
        end
        return (rt, out, "")
    catch e
        return (Any, out, string(e))
    end
end

function report(tag, f, tt)
    rt, bad, err = scan_unstable(f, tt)
    println("\n", "#"^70)
    println("# ", tag)
    println("#   sig = ", f, "  ", tt)
    println("#   inferred return type : ", rt, isconcretetype(rt) ? "  [concrete]" : "  [<-- NON-CONCRETE]")
    if !isempty(err); println("#   note: ", err); end
    println("#   non-concrete SSA values: ", length(bad))
    # show up to 12 distinct abstract types
    let seen = Set{Any}(), n = 0
        for (_, t) in bad
            t in seen && continue
            push!(seen, t); n += 1
            println("#     - ", t)
            n >= 12 && break
        end
    end
    println("#"^70)
    return warntype_str(f, tt)
end

KAT = typeof(ka)

# small probe: the KA field-access pattern at the top of every flux!
getflux(ka) = ka.kinfo.config.solver.flux
getkinfo(ka) = ka.kinfo
getfield_kdata(ka) = ka.kdata.field.faces

io = open("diag_warntype_rank$RANK.txt", "w")
function dump(title, s)
    println(io, "\n", "="^78, "\n== ", title, "\n", "="^78)
    println(io, s)
end

dump("PROBE getkinfo(ka)  [ka.kinfo field access]", report("getkinfo(ka)::ka.kinfo", getkinfo, (KAT,)))
dump("PROBE getflux(ka)   [ka.kinfo.config.solver.flux]", report("getflux(ka)", getflux, (KAT,)))
dump("PROBE getfield_kdata(ka) [ka.kdata.field.faces]", report("getfield_kdata(ka)", getfield_kdata, (KAT,)))

# ---- top-level flux! for the 4 face types (type-based; valid w/o instances) ----
for FT in (KitAMR.DomainFace, KitAMR.FullFace, KitAMR.HangingFace, KitAMR.BackHangingFace)
    # concrete parametrization present in this 2D/2F run
    CT = try FT{2,2} catch; FT{2,2,KitAMR.Maxwellian} end
    tag = "flux!($(nameof(FT)){...}, ka)"
    s = report(tag, KitAMR.flux!, (CT, KAT))
    dump("TOPLEVEL " * tag, s)
end

# ---- descend into inner calls using REAL objects (Cthulhu-style) ----
ff_i = findfirst(f -> isa(f, KitAMR.FullFace), faces)
df_i = findfirst(f -> isa(f, KitAMR.DomainFace), faces)
hf_i = findfirst(f -> isa(f, KitAMR.HangingFace), faces)

if ff_i !== nothing
    ff = faces[ff_i]
    s1 = report("make_face_vs(FullFace)", KitAMR.make_face_vs, (typeof(ff),))
    dump("INNER make_face_vs(FullFace)", s1)
    here_vs, there_vs = KitAMR.make_face_vs(ff)
    println("\n--- typeof(here_vs) = ", typeof(here_vs))
    s2 = report("calc_flux(CAIDVM, here_vs, there_vs, FullFace, ka)",
        KitAMR.calc_flux, (Type{CAIDVM}, typeof(here_vs), typeof(there_vs), typeof(ff), KAT))
    dump("INNER calc_flux(CAIDVM, FullFace)", s2)
    # update_flux!
    flux, micro = KitAMR.calc_flux(CAIDVM, here_vs, there_vs, ff, ka)
    s3 = report("update_flux!(flux,micro,FullFace,heavi)",
        KitAMR.update_flux!, (typeof(flux), typeof(micro), typeof(ff), typeof(here_vs.heavi)))
    dump("INNER update_flux!(FullFace)", s3)
end

if df_i !== nothing
    df = faces[df_i]
    s1 = report("make_face_vs(DomainFace)", KitAMR.make_face_vs, (typeof(df),))
    dump("INNER make_face_vs(DomainFace)", s1)
    here_vs = KitAMR.make_face_vs(df)
    s2 = report("calc_domain_flux(CAIDVM, here_vs, DomainFace, ka)",
        KitAMR.calc_domain_flux, (Type{CAIDVM}, typeof(here_vs), typeof(df), KAT))
    dump("INNER calc_domain_flux(CAIDVM, DomainFace)", s2)
end

if hf_i !== nothing
    hf = faces[hf_i]
    s1 = report("make_face_vs(HangingFace)", KitAMR.make_face_vs, (typeof(hf),))
    dump("INNER make_face_vs(HangingFace)", s1)
end

close(io)
println("\nWROTE diag_warntype_rank$RANK.txt")
finalize!(p4est, ka)
MPI.Finalize()
