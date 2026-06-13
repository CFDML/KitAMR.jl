# Output configuration
KitAMR.jl outputs the field as `.jld2` and `.pvtu` files for post-processing.
```@docs
Output
```
Supported vtk cell types are
```@docs
Pixel
```
```@docs
Voxel
```
```@docs
Triangle
```
```@docs
Tetra
```
Note that the vertices in `.pvtu` file are redundant. To obtain contour of the field, a cleaning process in post-process software is required. For example, using [Paraview](https://www.paraview.org), one needs to execute `Clean to Grid` before contour.

`Triangle` and `Tetra` are chosen for smoother field visualization, but will not improve the accuracy of the solution. The variables on their vertices are linearly reconstructed. 

