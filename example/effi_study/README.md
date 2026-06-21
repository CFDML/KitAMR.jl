# R2D Vacuum-Jet AMR Efficiency Study

This directory contains a cleaned benchmark derived from
`dev/example/r2d/r2d_vacuum_jet.jl`.

The benchmark keeps the active default vacuum-jet case:

- left aperture jet with open exterior,
- no bottom jet,
- no fixed-vacuum left-boundary branch,
- no initial packet branch,
- no velocity-space probe output during timing.

Four modes are submitted as independent jobs:

- `both_on.jl`: physical-space AMR + velocity-space AMR.
- `ps_only.jl`: physical-space AMR + uniform finest velocity grid.
- `vs_only.jl`: uniform finest physical grid + velocity-space AMR.
- `both_off.jl`: uniform finest physical and velocity grids.

The Slurm wrappers resolve the project root from their own location, so they are
portable after copying the repository to a cluster. Keep `example/effi_study`
and `slurm_templates` in the same project tree.

Defaults:

- `VJE_MAX_SIM_TIME=0.2`
- `VJE_ANIM_DT=0.05`
- `VJE_PS_BASE=16`
- `VJE_PS_LMAX=3`
- `VJE_VS_BASE=8`
- `VJE_VS_LMAX=4`

Submit all four modes separately:

```text
SBATCH_EXTRA="-p <partition> -A <account>" JULIA_SYSIMAGE=/path/to/sys.so SRUN_MPI=pmix bash example/effi_study/slurm/submit_modes.sh
```

Run analysis after the four jobs finish:

```text
JULIA_SYSIMAGE=/path/to/sys.so bash example/effi_study/slurm/submit_analysis.sh
```

Outputs are written under `example/effi_study/out` by default.  Analysis
writes `amr_effi_report.csv` and `summary.tsv`.
