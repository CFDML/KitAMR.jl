# KitAMR.jl

## Motivation
**KitAMR.jl** is a distributed adaptive Cartesian grid solver for kinetic equations, developed based on [P4est.jl](https://github.com/trixi-framework/P4est.jl). Its goal is to achieve large-scale parallel solving of 2D and 3D flows across all regimes, leveraging GPUs, differentiable programming, and machine learning to enhance solving efficiency.

Currently, the development is in its initial stages, with the following functionalities implemented:
- 2D flow
- Square domain grid generation
- Geometric adaptation of physical space grids
- Dynamic adaptation of physical space grids
- Dynamic adaptation of velocity space grids
- Load balancing based on physical space grid density
- VanLeer reconstruction
- UGKS fluxes
- Post-processing

## Theory
### Kinetic Methods
In kinetic theory, methods simulate the macroscopic motion of fluids by describing the evolution of distribution functions of fluid particles in velocity space over time. Such methods include direct simulation Monte Carlo (DSMC), discrete velocity method (DVM), lattice Boltzmann equation (LBE), gas-kinetic scheme, semi-Lagrangian method, implicit-explicit (IMEX) method, and others.
### Distribution Function
In kinetic methods, the state of the fluid is described by a distribution function $f(\mathbf{x},\mathbf{u},t)$. The distribution function represents the number of particles at a given moment in time at a specific physical space point in a particular velocity space element. Its normalization property is

$$
\rho(\mathbf{x},t)=\iiint_{-\infty}^{\infty}f(\mathbf{x},\mathbf{u},t)d\mathbf{u}
$$

, where $\rho(\mathbf{x},t)$ is the density of flows at $\mathbf{x}$ in time $t$. At higher temperatures, molecular internal degrees of freedom are excited, and the distribution function will depend on these degrees of freedom as well. For simplicity, we do not currently consider this situation.

<div align="center">
    <img src="https://i.postimg.cc/pTsjgW3T/20240703151108.jpg)](https://postimg.cc/ygSWRH14" alt="Image 1" width="300" title="distribution function of 2-D" style="margin-right: 1px;">
</div>
With the distribution function, macroscopic conservative variables can be caculated as follows:

$$
\begin{split}
&\rho = \iiint_{-\infty}^{\infty}f(\mathbf{x},\mathbf{u},t)d\mathbf{u}\\
&\rho\mathbf U = \iiint_{-\infty}^{\infty}\mathbf uf(\mathbf{x},\mathbf{u},t)d\mathbf{u}\\
&\rho E = \iiint_{-\infty}^{\infty} \frac 12 \mathbf u^2f(\mathbf{x},\mathbf{u},t)d\mathbf{u}\\
&\mathbf q=\iiint_{-\infty}^{\infty} \mathbf u\frac 12 \mathbf u^2f(\mathbf{x},\mathbf{u},t)d\mathbf{u}
\end{split}
$$

### Boltzmann Equation
The evolution of the distribution function follows the Boltzmann equation.

$$
\frac{\partial f}{\partial t}+\mathbf{u}\cdot\frac{\partial f}{\partial \mathbf{x}}+\mathbf{F}\cdot\frac{\partial f}{\partial \mathbf{u}}=\iiint_{-\infty}^{\infty}\int_0^{4\pi}(f^{\*}f_1^{\*}-ff_1)c_r\sigma d\Omega dc_1
$$

<!-- The Boltzmann equation holds a central position in kinetic theory. In the free molecular flow regime, the collisionless Boltzmann equation is applied, or its equilibrium solution, the Maxwellian distribution, is directly used. In the slip flow regime, the Boltzmann equation is expanded via the Chapman-Enskog method to obtain the first-order approximation, the Navier-Stokes equations, or the second-order approximation, the Burnett equations. Alternatively, a systematic asymptotic expansion yields a set of fluid mechanics equations. In the transition regime, solving flow problems involves approaches derived from the Boltzmann equation or its equivalent methods. -->
The Boltzmann equation holds a central position in kinetic theory, which is an integro-differential equation, with a four-fold integral term on the right-hand side. The complexity introduced by this term makes solving the Boltzmann equation challenging. Therefore, as a simplification, kinetic methods often solve its model equations instead.
### Shakhov Model Equation
The simplified equation obtained after simplifying the collision term is known as the Shakhov model equation:

$$
\begin{cases}
\frac{\partial f}{\partial t}+\mathbf{u}\cdot\frac{\partial f}{\partial \mathbf{x}}+\mathbf{F}\cdot\frac{\partial f}{\partial \mathbf{u}}=\frac{g^+-f}{\tau}\\
g^+=g[1+\frac45(1-\mathrm{Pr})\lambda^2\frac{\mathbf{u}\cdot\mathbf{q}}{\rho}(2\lambda \mathbf{u}^2-5)]
\end{cases}
$$

, Where $\mathrm {Pr}$ is the Prandtl number, $\mathbf{q}$ is the heat flux, and $\tau$ is the relaxation time typically given through phenomenological models to match experimental data. The Shakhov model can provide results close to the Boltzmann equation for mass, momentum, energy, and heat fluxes in situations where flow velocities are not very high and deviations from equilibrium are not very large. And the current solver is based on Shakhov model.

## Computing Methods
### Finite Volume Method (FVM)
The underlying discretization of the Boltzmann equation is based on finite volume method, with its discrete form as follows:

$$
(\bar{U}_j\Omega_j)^{n+1} = (\bar{U}_j \Omega_j)^n-\Delta t\sum_{\partial\Omega}\mathbf{F}^{\*}\cdot\Delta\mathbf{S}+\Delta t \bar{Q}_j\Omega_j
$$

, where $\bar U_j^n\coloneqq\frac{1}{\Omega_j}\int_{\Omega_j}Ud\Omega_j|^n$, $\bar Q_j^n\coloneqq\frac{1}{\Omega_j}\int_{\Omega_j}Qd\Omega_j$ are the cell-averaged conservative variable and source. And $\bar Q_j^{\*}$ and $\mathbf{F}^{\*}$ are respectively cell- and time-averaged sources and numerical flux, which are defined as $\bar Q_j^{\*}\coloneqq\frac{1}{\Delta t}\int_n^{n+1}\bar Q_jdt$ and $\mathbf F^{\*}\cdot \Delta \mathbf{S}\coloneqq \frac{1}{\Delta t}\int_n^{n+1}\mathbf F\cdot\Delta \mathbf{S}dt$.
It is an exact relation for the time evolution of the space averaged conservative avriables over cell $j$ from time step $n$ to $n+1$. And the numerical flux $\mathbf F^{\*}$ identifies completely a scheme by the way it approximates the time-averaged physical flux along each cell face.
### Numerical Flux
The numerical scheme that the solver currently supports is UGKS, with interface variables obtained from VanLeer reconstruction. For specific implementation details, please refer to the relevant literature.

### Discrete Velocity Method
DVM is one of the popular approaches for solving rarefied flow problems. In this method, the particle velocity space is discretized into a finite set of points and the numerical quadrature rule is utilized to approximate the integration of moments. Since the particle velocity space is discretized into a finite set of points, the continuum Boltzmann equation is reduced to the corresponding discrete velocity Boltzmann equation (DVBE). A great variety of algorithms have been developed, including unified gas-kinetic scheme (UGKS), discrete unified gas-kinetic scheme (DUGKS), semi-Lagrangian method, etc. Overall, most of the above methods can be applied from free molecular regime to continuum regime. But in order to make the quadrature error to be small enough, a large number of discrete velocity points are usually required. In particular, for fluid flows near continuum regime, the computational cost of DVM is much larger than those traditional CFD methods based on the Navier-Stokes equation.

### Adaptive Mesh Refinement (AMR)
Considering the above limitations of DVM, we adopt Adaptive Mesh Refinement (AMR) to improve solving efficiency. AMR reallocates computational resources based on flow features, balancing efficiency and accuracy. 

<div align="center">
    <img src="https://i.postimg.cc/FzQx9cSD/PV-AMR.png" alt="Image 1" title="Simultaneous AMR" width="700" style="margin-right: 1px;">
</div>

<div align="center">
    <img src="https://i.postimg.cc/7Lvp4xWX/cegur-0qkdq.gif" alt="Image 1" title="physical space" width="300" style="margin-right: 1px;">
    <img src="https://i.postimg.cc/rwVV31C8/d5ubx-sf3vl.gif" alt="Image 2" title="simultaneously, in velocity space" width = "300"style="margin-left: 1px;">
</div>
<!-- Generally speaking, the criteria for refining or coarsening physical space grids are based on the magnitude of macroscopic variable gradients. This approach aims to use finer grids to capture more pronounced variations, thereby modeling more complex flow fields. On the other hand, criteria for velocity space grids are based on the proportion of energy of particles within the grid relative to the total energy of all particles at that physical point. This is because particles with higher energy proportions exert a more significant influence on the macroscopic behavior of the flow, necessitating more accurate numerical integration.  -->
The current solver discretizes with tree-based Cartesian grids, adapting the mesh in both physical and velocity space simultaneously, according to characteristics including velocity field gradients, boundary complexity, and energy proportion at phase points. Users can also customize the criteria for grid refinement or coarsening according to their own needs.

### Cartesian Grids
Cartesian grids offer advantages such as high grid quality and automation but struggle with complex boundary shapes. In the subsequent development of the solver, we will introduce Immersed Boundary Method (IBM) to address this challenge.

## Installation
The current solver does not yet provide a user interface. If you'd like to try it out, we recommend cloning the code repository from Git to your local machine and then using the Julia package manager to install its dependencies in the directory where the code resides:
```julia
julia> ]
(v1.10) pkg> activate .
(//your_path//) pkg> instantiate
```
Soon, we will provide a complete user interface with customizable boundaries, and we will register it in the [Julia package registry](https://github.com/JuliaRegistries/General). Stay tuned!
## Current Usage
Currently, flow properties, initial contitions and boundary conditions are defined in `Global_Data` in `/src/types.jl`. It should be noted specifically that the boundary follows the sequence of -x, +x, -y, +y. The maximum level of the refinement is defined in `/src/abstract.jl` by `DVM_PS_MAXLEVEL` and `DVM_VS_MAXLEVEL` respectively. Grid refinement and coarsening behaviors are defined in `adaptive.jl` and `vs_adaptive.jl`. Feel free to modify them according to your needs.