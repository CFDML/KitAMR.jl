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

function check!(i,ps4est,amr)
    if amr.global_data.status.residual.step%RES_CHECK_INTERVAL==0
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            i+=1
            @show i
            sim_time = amr.global_data.status.sim_time
            @show sim_time
            res = maximum(amr.global_data.status.residual.residual)
            @show res
            pp = PointerWrapper(ps4est)
            local_num_quadrants = pp.local_num_quadrants[]
            @show local_num_quadrants
            ref_vs_num = amr.global_data.status.max_vs_num
            @show ref_vs_num
        end
    end
    check_for_save!(ps4est,amr)
end