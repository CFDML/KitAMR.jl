# Cluster deployment (custom MPI & p4est)

By default KitAMR needs **no** special setup: `] add KitAMR` pulls in
[`P4est.jl`](https://github.com/trixi-framework/P4est.jl) with the prebuilt
`P4est_jll` binary and `MPICH_jll`, which work out of the box on laptops and
workstations.

On an HPC cluster you usually want the cluster's **own MPI** (InfiniBand /
Slurm-integrated, etc.). Because the prebuilt `P4est_jll` is compiled against
`MPICH_jll` — not your cluster MPI — and P4est.jl does **not** support mixing a
system MPI with a JLL `p4est`, you must use:

* the cluster's system MPI, **and**
* a `p4est` built against that *same* MPI.

This is a **one-time, per-environment** setup (not per run) and is configured
entirely through [`Preferences.jl`](https://github.com/JuliaPackaging/Preferences.jl)
and [`MPIPreferences.jl`](https://juliaparallel.org/MPI.jl/stable/configuration/) —
the same mechanism documented by MPI.jl and P4est.jl. Nothing is patched inside
the installed package; everything lands in your project's `LocalPreferences.toml`.

## Step 1 — Build `p4est` against the cluster MPI

Download a `p4est` source release whose version matches the JLL major/minor. Run the following script in the p4est source tree:

```bash
module load <your-mpi>                 # provides mpicc (and libmpi at run time)
PREFIX="$PWD/../libp4est"              # absolute install prefix
./bootstrap
./configure CC=mpicc --enable-mpi --enable-shared \
    --prefix="$PREFIX" LDFLAGS="-Wl,-rpath,$PREFIX/lib"
make -j"$(nproc)"
make install
```

The `-Wl,-rpath,$PREFIX/lib` flag embeds the library directory into
`libp4est.so` so it can locate `libsc.so.3` at run time **without**
`LD_LIBRARY_PATH`. After `make install` you should have
`$PREFIX/lib/{libp4est.so, libsc.so, libsc.so.3, …}`.

!!! tip
    If your cluster already provides a system-MPI `p4est` (e.g. `module load
    p4est`), skip this step and use that module's `lib/` directory in Step 2.

## Step 2 — Point P4est at it (Preferences)

Run this in the project that depends on KitAMR. **Do not** `using KitAMR` /
`using P4est` here — that would try to load the very library you are about to
configure (and fail if the current backend is broken). Only the lightweight
`Preferences` / `MPIPreferences` packages are loaded:

```julia
using Preferences, UUIDs, MPIPreferences
const P4EST = UUID("7d669430-f675-4ae7-b43e-fab78ec5a902")   # registered P4est.jl (stable uuid)
lib = "/abs/path/to/libp4est/lib"                            # The path to libp4est is the same as the PREFIX in Step 1 (or the p4est module)
set_preferences!(P4EST,
    "libp4est" => joinpath(lib, "libp4est.so"),
    "libsc"    => joinpath(lib, "libsc.so");
    force = true)
MPIPreferences.use_system_binary()
```

This writes your project's `LocalPreferences.toml` — a per-project, per-machine
file that is **not** part of any package. The p4est preference is keyed under
P4est's UUID (this is exactly the upstream P4est.jl procedure); the MPI choice is
handled by MPIPreferences.

## Step 3 — Rebuild MPI and restart

```bash
julia --project -e 'using Pkg; Pkg.build("MPI")'
```

Then **restart Julia**. These are load-time preferences, so P4est/KitAMR are
recompiled and pick up the new backend on the next `using KitAMR`.

## Running jobs

With the rpath from Step 1, `p4est`/`libsc` need no extra environment — just make
the MPI runtime available:

```bash
module load <your-mpi>
mpiexec -n <N> julia --project your_script.jl
```

If you did **not** build with rpath (or use a prebuilt p4est that lacks one),
export the library directory **before** launching Julia (setting it inside a
running REPL is too late):

```bash
export LD_LIBRARY_PATH=/abs/path/to/libp4est/lib:$LD_LIBRARY_PATH
```

## Faster startup & lower memory: build a sysimage

On a cluster you launch **one Julia process per MPI rank**, and every rank pays the
full Julia runtime cost on its own: it JIT-compiles the same methods and keeps a
*private* copy of that native code, the compiler's data structures and the GC heap.
This private part is roughly **two-thirds (~70%) of each rank's baseline memory**,
and it is replicated on every rank — at 64 ranks/node that is a large fixed cost
before a single phase-space cell is stored, plus a multi-second JIT delay at every
launch.

A **sysimage** — a `.so` produced by
[`PackageCompiler.jl`](https://github.com/JuliaLang/PackageCompiler.jl) — bakes
KitAMR's compiled code into a memory-mapped image. The OS maps that image **once per
node and shares it across all ranks**, so the code moves from private to shared
memory. In a controlled measurement this cut each rank's private memory by **~50%**
and roughly **halved the true per-node Julia baseline**, and it removes the startup
JIT. The parallel (pure-MPI) structure is unchanged — this is purely a launch-time
image swap.

### Step 1 — Install PackageCompiler (a build tool)

It is only needed to *build* the image, not to run KitAMR. Install it in your
default environment so it stays separate from the environment you run jobs in:

```bash
julia -e 'using Pkg; Pkg.add("PackageCompiler")'
```

With the default `JULIA_LOAD_PATH`, a later `julia --project` still finds
PackageCompiler through the default environment while KitAMR resolves from your
project environment.

### Step 2 — Write a precompile workload

A small script that exercises the code paths you actually run (initialize → AMR →
flux/iterate → output). PackageCompiler records every method specialization it
touches and bakes them in. Keep it small and gentle — what matters is the *types and
code paths*, not the physics or the problem size (see the **Coverage** section
below):

```julia
# precompile_workload.jl
using KitAMR, MPI
MPI.Init()
ic(mp, kinfo) = [1.0, 0.0, 0.0, 1.0]                 # any valid IC; only the types matter
solver = Solver(; DIM=2, NDF=2, CFL=0.4, AMR_PS_MAXLEVEL=1, AMR_VS_MAXLEVEL=3,
    PS_DYNAMIC_AMR=true, VS_DYNAMIC_AMR=true, flux=CAIDVM,
    time_marching=CIP_Marching, max_sim_time=1.0)
gas = Gas(; K=1.0, Kn=1e-3, ω=0.81, ωᵣ=0.81)
config = Configure(solver; geometry=[0.,1.,0.,1.], trees_num=[8,8],
    quadrature=[-10.,10.,-10.,10.], vs_trees_num=[8,8], IC=PCoordFn(ic),
    domain=[Domain(Period,i) for i in 1:4],
    output=Output(solver), gas=gas, user_defined=UDF())
p4est, ka = initialize(config; prerefine_steps=1)
solve!(p4est, ka; max_steps=3, listen_for_save=false, progress=false)
save_result(p4est, ka; dir_path=joinpath(mktempdir(), "r"))
KitAMR.finalize!(p4est, ka)
MPI.Finalize()
```

### Step 3 — Build the sysimage

```julia
# build_sysimage.jl
using PackageCompiler
create_sysimage(
    [:KitAMR, :MPI, :JLD2];
    sysimage_path = "kitamr_sys.so",
    precompile_execution_file  = "precompile_workload.jl",
    # precompile_statements_file = "stmts.jl",       # optional, recommended — see Coverage
    sysimage_build_args = `--strip-metadata`,        # smaller image, no debug info
    # cpu_target = "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)",  # for cross-CPU portability
)
```

```bash
julia --project build_sysimage.jl                    # takes several minutes -> kitamr_sys.so
```

!!! warning
    Build on the **same CPU microarchitecture** you run on (the baked native code is
    CPU-specific), or set `cpu_target`. **Rebuild** whenever you update KitAMR or its
    dependencies — a stale image silently runs old code.

### Step 4 — Run with the sysimage

Point Julia at it; nothing else changes:

```bash
mpiexec -n <N> julia --project --sysimage=kitamr_sys.so your_script.jl
# SLURM:
srun    -n <N> julia --project --sysimage=kitamr_sys.so your_script.jl
```

!!! tip
    Pair it with `--heap-size-hint=<size>` (e.g. `--heap-size-hint=3G`) to cap the GC
    heap and keep peak RSS lower — the cheapest guard against transient
    out-of-memory on memory-bound runs.

### Coverage — what gets baked in, and what if a method is missed

A sysimage is an **optimization, not a sealed set of methods**. Any method *not*
captured by the workload is simply JIT-compiled the normal way the first time it is
called at run time. Results stay correct; you only lose the benefit *for that
method* — its code stays private per rank, and there is a one-time compilation pause
at first use. So **100% coverage is not required**, and a miss never causes a crash.

What a method gets compiled for is the **types and code paths** exercised, **not the
problem size**. Larger meshes, more ranks (beyond two), deeper AMR levels or more
steps are run-time *data* — they add no new compiled methods. To cover what a
production run needs, the workload must match the same *type combinations* and
*trigger* the same paths, at any small scale:

* **Types** — `DIM`/`NDF`, `flux`/`time_marching`, the boundary / immersed-boundary
  kinds, and the output cell types. A 2-D image does **not** cover 3-D; a periodic
  workload does **not** cover immersed-boundary methods.
* **Paths** — dynamic AMR (make it both **refine and coarsen**), and the **parallel
  communication paths** (ghost exchange, partition). The latter are only invoked
  with **≥ 2 ranks**, so a single-process build misses them.

Because `precompile_execution_file` runs single-process, capture the parallel
methods by also recording a `precompile_statements_file` from a small **2-rank** run
of a representative case, then pass it to `create_sysimage` (Step 3):

```bash
mpiexec -n 2 julia --project --trace-compile=stmts.jl your_small_case.jl
```

To check for residual gaps, run a small case **with** the finished sysimage and
`--trace-compile`; anything still printed is a method that was *not* in the image —
feed it back and rebuild:

```bash
mpiexec -n 2 julia --project --sysimage=kitamr_sys.so --trace-compile=missing.jl your_small_case.jl
```

[`SnoopCompile.jl`](https://github.com/timholy/SnoopCompile.jl) automates this
discovery if you want a more systematic sweep.

## Switching back to the bundled binaries

```julia
using Preferences, UUIDs, MPIPreferences
p4est_uuid = UUID("7d669430-f675-4ae7-b43e-fab78ec5a902")
delete_preferences!(p4est_uuid, "libp4est", "libsc"; force = true)
MPIPreferences.use_jll_binary()
```

Then `Pkg.build("MPI")` and restart Julia.

## Troubleshooting

* **`could not load library ".../libp4est.so" … libsc.so.3: cannot open shared
  object file: No such file or directory`** — the path preference is fine
  (`libp4est.so` was found), but the dynamic linker cannot resolve its dependency
  `libsc.so.3`. Fix with any of:
    * rebuild with `-Wl,-rpath,<lib>` (Step 1), or
    * `patchelf --set-rpath <lib> <lib>/libp4est.so`, or
    * `export LD_LIBRARY_PATH=<lib>:$LD_LIBRARY_PATH` before starting Julia.

    Diagnose with `ldd <lib>/libp4est.so` (look for `libsc.so.3 => not found`).

* **Warning that a system MPI is used with a JLL `p4est` (or vice versa)** — the
  two must be set together. Make sure you ran *both* `set_preferences!(P4EST, …)`
  and `MPIPreferences.use_system_binary()`.

* **`Detected version … of p4est … we only support v2.x from v2.3.0 on`** — build
  a 2.x `p4est` (match `P4est_jll` 2.8.x).
