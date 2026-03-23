using MPI
MPI.Init()
mutable struct MeshData
    in_box::Int
    in_solid::Bool
    in_search_radius::Int # In i-th ib's search radius. If not, i is set to 0.
    is_ghost_cell::Bool
end
function MeshData()
    return MeshData(0, false, 0, false)
end

function create_MPI_static_data(::Type{T}) where {T}
    types = fieldtypes(T)
    blocklengths = ones(Cint, length(types))
    displacements = [Cint(fieldoffset(T, i)) for i = 1:length(types)]
    datatypes = [MPI.Datatype(t) for t in types]
    mpi_type = MPI.Types.create_struct(blocklengths, displacements, datatypes)
    MPI.Types.commit!(mpi_type)
end
MPI_type = create_MPI_static_data(MeshData)
MPI.Types.size(MPI_type)

if MPI.Comm_rank(MPI.COMM_WORLD)==0
    mesh_datas = [MeshData(), MeshData(), MeshData()]
    buffer = MPI.Buffer(mesh_datas, 2, MPI_type)
    req = MPI.Isend(buffer, MPI.COMM_WORLD; dest = 1, tag = 630)
    MPI.Wait(req)
else
    # rec_vec = Vector{MeshData}(undef,2)
    rec_buffer = MPI.Buffer(MeshData(), 1, MPI_type)
    req = MPI.Irecv!(rec_buffer, MPI.COMM_WORLD; source = 0, tag = 630)
    MPI.Wait(req)
    @show rec_buffer
    # mesh_data = rec_vec[1]
    # @show mesh_data.is_ghost_cell
end
using StructArrays
mesh_datas = [MeshData(), MeshData(), MeshData()]
datas = StructArray(mesh_datas)
datas.in_solid
datas.in_solid[1] = true
mesh_datas[1]==datas[1]
```
From the document of StructArrays.jl:
    Finally, when created from an array-of-structs, StructArrays creates a copy of the "parent" data. This effectively "detaches" 
    the StructArray from the original data.
```
