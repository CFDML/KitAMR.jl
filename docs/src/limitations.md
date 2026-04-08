# Limitations

In current stage, KitAMR.jl has several noticeable limitations, including
- Flows with volume forces (which will introduce derivatives with respect to velocity coordinates) are not supported.
- Only translational non-equilibrium effects is introduced.
- Only support fully diffuse Maxwellian gas-surface interaction model.
- Restart is not supported in consideration of the enormous size of data in complete phase space.
- Only support Euler method for time marching, which only exhibits first order accuracy. For transient problems, higher order time marching scheme is required to be consistent with spatial and quadrature accuracy.

Any contribution is welcomed. Please feel free to ask us questions and chat with us at any time if you're unsure about anything.