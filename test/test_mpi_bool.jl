using MPI
MPI.Init()
COMM_DATA_TAG = 10000
N = MPI.Comm_size(MPI.COMM_WORLD)
buffer = [Ref(false) for _ = 1:N]
reqs = Vector{MPI.Request}(undef, 0)
rank = MPI.Comm_rank(MPI.COMM_WORLD)
buffer[rank+1][] = isodd(rank)
for i = 1:N
    i-1==rank&&continue
    sreq = MPI.Isend(buffer[rank+1], MPI.COMM_WORLD; dest = i-1, tag = COMM_DATA_TAG+rank)
    push!(reqs, sreq)
end
for i = 1:N
    i-1==rank&&continue
    rreq = MPI.Irecv!(buffer[i], MPI.COMM_WORLD; source = i-1, tag = COMM_DATA_TAG+i-1)
    push!(reqs, rreq)
end
MPI.Waitall(reqs)
flags = [x[] for x in buffer]
if rank==0
    @show flags
end
MPI.Finalize()
