# Parallel

During a simulation with distributed parallel, the communication among processors is crucial. In KitAMR.jl, the communication is supported by [MPI.jl](https://github.com/JuliaParallel/MPI.jl).

## Ghost layers

KitAMR.jl decomposes the computational domain in physical space. With explicit spatial discretization, the update of a cell only depends on its neighboring cells. If some of these neighboring cells lie on other processors, they are locally marked as ghost. The variables of the ghost cells are updated by MPI communication. Related functions are
```@docs
data_exchange!
```
```@docs
slope_exchange!
```
```@docs
solid_exchange!
```
In contrast to traditional CFD solver, KitAMR.jl also discretizes velocity space. The resolution of this discretization varies as time and physical coordinates changes. This results in the size of communicated data changing. Currently, the size of the communicated data is unified to hold the largest velocity space that exists in ghost layers. The size is obtained by
```@docs
get_vs_num
```
Whether a unified communication size decreases the efficiency noticeably still requires testing.

## Partition

After an AMR process, the grids density on different processors may changes a lot. To balance the load, partition can be performed by calling
```@docs
partition!
```
The grids in physical space are encoded as a 1-dimensional sequence by Morton code, and then partitioned according to the weight provided by
```@docs
partition_weight
```
The backend functions are provided by p4est, insuring high efficiency. 

When the partition is finished, a data exchange process is performed to transfer all of the data to where it is required. The communication only occurs to those partitioned cells.

All of the functionalities are wrapped into
```@docs
ps_partition!
```
