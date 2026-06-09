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
