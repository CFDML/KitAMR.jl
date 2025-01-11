function unpack(t)
    return (getfield(t,i) for i in 1:nfields(t))
end
function listen_for_save!(;rank=0)
    if MPI.Comm_rank(MPI.COMM_WORLD)==rank
        Base.start_reading(stdin)
    end
    return nothing
end
function check_for_save!(ps4est::P_pxest_t,amr;rank=0)
    save_flag = amr.global_data.status.save_flag
    if MPI.Comm_rank(MPI.COMM_WORLD)==rank
        if bytesavailable(stdin)!=0
            input = readline(stdin)
            if input=="save"
                save_flag[] = true
            end
        end
    end
    MPI.Bcast!(save_flag,rank,MPI.COMM_WORLD)
    if save_flag[] == true
        MPI.Comm_rank(MPI.COMM_WORLD)==rank&&println("saving...")
        save_result(ps4est,amr)
        MPI.Comm_rank(MPI.COMM_WORLD)==rank&&println("done.")
        save_flag[] = false
    end
    return nothing
end

function MPI_tune(f::Function,args...)
    for i in 1:MPI.Comm_size(MPI.COMM_WORLD)
        if MPI.Comm_rank(MPI.COMM_WORLD)==i-1
            f(args...)
        end
        MPI.Barrier(MPI.COMM_WORLD)
    end
end

