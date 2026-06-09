# Flux

Numerical fluxes in KitAMR.jl are abstracted as
```@docs
AbstractFluxType
```
To distinguish with Navier-Stokes and other types of fluxes (which may be implemented in future version), the fluxes corresponding to DVM are abstracted to
```@docs
AbstractDVMFluxType
```
Currently, the best supported flux is 
```@docs
CAIDVM
```

Here are functions related to flux calculation and update:

## Backend functions

```@docs
flux!(::KA)
```
---

```@docs
flux!(::Type{F},::DomainFace,::KA) where {F<:AbstractFluxType}
flux!(::Type{F},::FullFace,::KA) where {F<:AbstractFluxType}
flux!(::Type{F},::HangingFace{DIM,NDF},::KA) where {DIM,NDF,F<:AbstractFluxType}
flux!(::Type{F},::BackHangingFace{DIM,NDF},::KA) where {DIM,NDF,F<:AbstractFluxType}
```
Methods for different types of faces. The numerical-flux scheme `F` is passed as a type parameter (read once from the configuration by [`flux!(::KA)`](@ref) and threaded through the face loop) so that the scheme is resolved statically. Through these methods, different scenarios are decomposed into a unified interface, and the same core functions can be invoked for solving.

---
```@docs
update_flux!
```
A shell for the dispatch according to whether the macro flux update is required. The method for [`DomainFace`](@ref) is separately defined: [`update_domain_flux!`](@ref).

---
```@docs
update_domain_flux!
```
Refer to [`update_flux!`](@ref).

---
```@docs
update_micro_flux!
```
Update `flux` in [`AbstractVsData`](@ref). Mapping between velocity cells with different refinement level is considered.

---
```@docs
update_macro_flux!
```
Update `flux` in [`AbstractPsData`](@ref), which is necessary to predict the macroscopic state in implicit DVM. The conserved correction in the mapping between different velocity cells also depends on this.

---
```@docs
make_face_vs
```
Construct [`FaceVsData`](@ref) corresponding to `here_data` and `there_data` in [`AbstractFace`](@ref).

## User interfaces

Developers may hope to use their own numerical flux to meet their requirements. KitAMR.jl provides compact user interfaces. User only needs to define the flux across a single inner face and a single domain face. KitAMR.jl will automatically dispatch them for different mesh cases. As an example, the methods corresponding to [`CAIDVM`](@ref) are defined as
```@docs
calc_flux(::Type{CAIDVM},here_vs,there_vs,flux_data::Union{FullFace,FluxData},ka::KA{DIM,NDF}) where{DIM,NDF}
```
---
```@docs
calc_domain_flux(::Type{CAIDVM},::FaceVsData,::DomainFace{DIM,NDF,Maxwellian},::KA) where{DIM,NDF}
calc_domain_flux(::Type{CAIDVM},::FaceVsData,::DomainFace{DIM,NDF,SuperSonicInflow},::KA) where{DIM,NDF}
calc_domain_flux(::Type{CAIDVM},::FaceVsData,::DomainFace{DIM,NDF,UniformOutflow},::KA) where{DIM,NDF}
calc_domain_flux(::Type{CAIDVM},::FaceVsData,::DomainFace{2,NDF,InterpolatedOutflow},::KA) where{NDF}
```
Compute the flux across a single domain face. Methods for different types [`AbstractBoundCond`](@ref) of [`Domain`](@ref) are defined.

---

A positivity preserving reconstruction is adopted when compute [`CAIDVM`](@ref) flux:
```@docs
positivity_preserving_reconstruct
```