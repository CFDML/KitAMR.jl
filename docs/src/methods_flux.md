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
flux!(::KitAMR_Data)
```
---

```@docs
flux!(::DomainFace,::KitAMR_Data)
flux!(::FullFace,::KitAMR_Data)
flux!(::HangingFace{DIM,NDF},::KitAMR_Data) where{DIM,NDF}
flux!(::BackHangingFace{DIM,NDF},::KitAMR_Data) where{DIM,NDF}
```
Methods for different types of faces. Through these methods, different scenarios are decomposed into a unified interface, and the same core functions can be invoked for solving.

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
Construct [`Face_VS_Data`](@ref) corresponding to `here_data` and `there_data` in [`AbstractFace`](@ref).

## User interfaces

Developers may hope to use their own numerical flux to meet their requirements. KitAMR.jl provides compact user interfaces. User only needs to define the flux across a single inner face and a single domain face. KitAMR.jl will automatically dispatch them for different mesh cases. As an example, the methods corresponding to [`CAIDVM`](@ref) are defined as
```@docs
calc_flux(::Type{CAIDVM},here_vs,there_vs,flux_data::Union{FullFace,FluxData},amr::KitAMR_Data{DIM,NDF}) where{DIM,NDF}
```
---
```@docs
calc_domain_flux(::Type{CAIDVM},::Face_VS_Data,::DomainFace{DIM,NDF,Maxwellian},::KitAMR_Data) where{DIM,NDF}
calc_domain_flux(::Type{CAIDVM},::Face_VS_Data,::DomainFace{DIM,NDF,SuperSonicInflow},::KitAMR_Data) where{DIM,NDF}
calc_domain_flux(::Type{CAIDVM},::Face_VS_Data,::DomainFace{DIM,NDF,UniformOutflow},::KitAMR_Data) where{DIM,NDF}
calc_domain_flux(::Type{CAIDVM},::Face_VS_Data,::DomainFace{2,NDF,InterpolatedOutflow},::KitAMR_Data) where{NDF}
```
Compute the flux across a single domain face. Methods for different types [`AbstractBoundCondType`](@ref) of [`Domain`](@ref) are defined.

---

A positivity preserving reconstruction is adopted when compute [`CAIDVM`](@ref) flux:
```@docs
positivity_preserving_reconstruct
```