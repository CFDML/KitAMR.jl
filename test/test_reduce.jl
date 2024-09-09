using MPI
MPI.Init()
buffer = MPI.Comm_rank(MPI.COMM_WORLD)*ones(2)
if MPI.Comm_rank(MPI.COMM_WORLD)==0
    buffer[1] = MPI.Reduce(
        buffer[1],
        (x,y)->(
            max(x,y)
        ),
        MPI.COMM_WORLD
    )
    buffer[2] = MPI.Reduce(
        buffer[2],
        (x,y)->(
            x+y
        ),
        MPI.COMM_WORLD
    )
    @show buffer
else
    MPI.Reduce(
        buffer[1],
        (x,y)->(
            max(x,y)
        ),
        MPI.COMM_WORLD
    )
    MPI.Reduce(
        buffer[2],
        (x,y)->(
            x+y
        ),
        MPI.COMM_WORLD
    )
end

MPI.Finalize()