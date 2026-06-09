# Limitations

In current stage, KitAMR.jl has several noticeable limitations, including
- Flows with volume forces (which will introduce derivatives with respect to velocity coordinates) are not supported.
- Only translational non-equilibrium effects is introduced.
- Only support fully diffuse Maxwellian gas-surface interaction model.
- Only support Euler method for time marching, which only exhibits first order accuracy. For transient problems, higher order time marching scheme is required to be consistent with spatial and quadrature accuracy.

Restart from a checkpoint is supported via [`save_for_restart`](@ref) / [`restart`](@ref), including
resuming on a different number of MPI ranks. Be aware that a checkpoint stores the full phase-space
state (every velocity cell's distribution function), so the files can be large — `save_for_restart`
reports the estimated size before writing.

Any contribution is welcomed. Please feel free to ask us questions and chat with us at any time if you're unsure about anything.