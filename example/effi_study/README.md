# R2D Vacuum-Jet AMR Efficiency Study

This directory contains a cleaned benchmark derived from
`dev/example/r2d/r2d_vacuum_jet.jl`.

The benchmark keeps the active default vacuum-jet case:

- left aperture jet with Maxwellian wall exterior at the vacuum-field temperature,
- no bottom jet,
- no fixed-vacuum left-boundary branch,
- no initial packet branch,
- no velocity-space probe output during timing.

The computational box is derived from the velocity range:

```text
L = VRANGE * 0.1
geometry = [0, L, -L, L]
```

The jet enters from the center of the left boundary, so the non-inlet boundaries are one fastest
travel distance away from the inlet point. The default physical roots are `8 x 16`, so root cells stay square in the
`1:2` box.  With `VJE_PS_LMAX=5`, the finest physical resolution is close to
the original `4 x 4`, `16 x 16`, level-5 setup.  Physical-space AMR uses a
custom UDF criterion that reuses the already-computed `ps_data.lohner` and
keeps only the temperature row, instead of the default density-plus-temperature
criterion.

Available entry points:

- `both_on.jl`: physical-space AMR + velocity-space AMR.
- `both_off.jl`: uniform finest physical and velocity grids.  The physical mesh
  is built from a lower root grid and statically refined to the finest
  resolution, so high-rank runs keep better p4est locality than a finest-as-root
  mesh.
- `both_on_3d.jl`: 3D physical-space AMR + velocity-space AMR.  The inlet is
  the same left-boundary jet, extended to a circular pipe on the `(y,z)` inlet
  plane.
- `both_off_3d.jl`: 3D uniform finest physical and velocity grids.  This is
  usually much larger than the 2D reference and is protected by the same dense
  estimate guard.

`both_off` is the dense reference and can be very large at the default finest
resolution. Set `VJE_ALLOW_DENSE=1` only after checking memory, or lower
`VJE_PS_LMAX` / `VJE_VS_LMAX` for local smoke tests.

For `both_off`, the root grid is chosen automatically from the MPI size: it is
the coarsest root grid whose tree count is at least
`VJE_BOTH_OFF_ROOTS_PER_RANK * np` and then uniformly refined by
`AMR_PS_MAXLEVEL` to recover the same finest physical resolution.  With the
default 2D reference on 1024 ranks this gives `32 x 64` roots and static level
3, instead of `256 x 512` roots.  Override the lower bound with
`VJE_BOTH_OFF_MIN_ROOTS` if needed.

The Slurm wrappers resolve the project root from their own location, so they are
portable after copying the repository to a cluster. Keep `example/effi_study`
and `slurm_templates` in the same project tree.

Defaults:

- `VJE_MAX_SIM_TIME=0.1`
- `VJE_ANIM_DT=0.025`
- `VJE_DOMAIN_TRAVEL_TIME=0.1`
- `VJE_PS_BASE_X=8`
- `VJE_PS_BASE_Y=16`
- `VJE_PS_BASE_Z=16` for `VJE_DIM=3`
- `VJE_PS_LMAX=5`
- `VJE_VS_BASE=8`
- `VJE_BOTH_OFF_ROOTS_PER_RANK=1`
- `VJE_BOTH_OFF_MIN_ROOTS=0`
- `VJE_VS_LMAX` is estimated from the lowest-temperature Maxwellian unless
  explicitly set. The default states give `VJE_VS_LMAX=5`, resolving one
  `sigma_min` with at least three finest velocity cells.

For `VJE_DIM=3`, the box is

```text
geometry = [0, L, -L, L, -L, L]
```

The default `8 x 16 x 16` roots keep `dx = dy = dz`, so the finest physical
resolution is the same as the 2D setup in every coordinate direction.

Submit the default 2D modes:

```text
SBATCH_EXTRA="-p <partition> -A <account>" JULIA_SYSIMAGE=/path/to/sys.so SRUN_MPI=pmix bash example/effi_study/slurm/submit_modes.sh
```

Submit the 3D AMR run:

```text
VJE_DIM=3 SBATCH_EXTRA="-p <partition> -A <account>" JULIA_SYSIMAGE=/path/to/sys.so SRUN_MPI=pmix bash example/effi_study/slurm/submit_modes.sh
```

By default, `VJE_DIM=3` submits only `both_on_3d`, because `both_off_3d` can be
very large at the matching finest resolution. To also submit the dense 3D
reference, set `VJE_MODES="both_on_3d both_off_3d"` and `VJE_ALLOW_DENSE=1`
after checking the memory estimate.

Run analysis after the selected jobs finish:

```text
JULIA_SYSIMAGE=/path/to/sys.so bash example/effi_study/slurm/submit_analysis.sh
```

Outputs are written under `example/effi_study/out` by default.  Analysis
writes `amr_effi_report.csv` and `summary.tsv`.
