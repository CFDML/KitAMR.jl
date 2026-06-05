"""
$(SIGNATURES)
Perform a check procedure. Simulation status will be output every `ST_CHECK_INTERVAL` steps.
"""
function check!(p4est,ka)
    if ka.kinfo.status.step%ka.kinfo.config.solver.ST_CHECK_INTERVAL==0
        execute_check!(p4est,ka)
        check_for_save!(p4est,ka)
    end
end
function execute_check!(p4est,ka)
    status = ka.kinfo.status
    max_vs_num,total_phase_num = check_vs_num(ka)
    status.max_vs_num = max_vs_num;status.total_phase_num = total_phase_num
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        println("Iteration: $(ka.kinfo.status.step)")
        sim_time = ka.kinfo.status.sim_time
        println("Simulation time: $sim_time")
        res = maximum(ka.kinfo.status.residual.residual)
        println("Residual: $res")
        println("Maximum number of velocity grids: $max_vs_num")
        pp = PointerWrapper(p4est)
        global_num_quadrants = pp.global_num_quadrants[]
        println("Total number of physical grids: $global_num_quadrants")
        println("Total number of phase grids: $total_phase_num")
    end
end
"""
$(SIGNATURES)
Start listening for the input from command line. Must be called before calling [`check_for_save!`](@ref).
"""
function listen_for_save!(;rank=0)
    if MPI.Comm_rank(MPI.COMM_WORLD)==rank
        Base.start_reading(stdin)
    end
    return nothing
end
"""
$(SIGNATURES)
Check whether `save` is input during the iteration. If the result is true, a saving process will be executed.
"""
function check_for_save!(p4est::P_pxest_t,ka;rank=0)
    save_flag = ka.kinfo.status.save_flag
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
        save_result(p4est,ka)
        MPI.Comm_rank(MPI.COMM_WORLD)==rank&&println("done.")
        save_flag[] = false
    end
    return nothing
end

function check_for_animsave!(p4est::P_pxest_t,ka;path="./animation")
    output = ka.kinfo.config.output
    output.anim_dt <= 0. && return nothing   # animation disabled; safe no-op
    sim_time = ka.kinfo.status.sim_time
    if sim_time==0.
        step = output.anim_index
        for (celltype,suffix) in celltype_outputs(output.vtk_celltype)
            save_anim_field(path,p4est,ka,celltype,step,sim_time,0,suffix)
        end
        save_anim_vs(path,ka,step,sim_time,0)
    elseif sim_time >= (output.anim_index+1)*output.anim_dt - 1e-10*output.anim_dt
        step = output.anim_index+1
        for (celltype,suffix) in celltype_outputs(output.vtk_celltype)
            save_anim_field(path,p4est,ka,celltype,step,sim_time,1,suffix)
        end
        save_anim_vs(path,ka,step,sim_time,1)
        output.anim_index += 1
    end
end

function check_vs_num(ka::KA)
    trees = ka.kdata.field.trees.data
    buffer = zeros(Int,2) # [max,total]
    for tree in trees
        for ps_data in tree
            isa(ps_data,InsideSolidData)&&continue
            vs_num = ps_data.vs_data.vs_num
            buffer[1] = max(buffer[1],vs_num)
            buffer[2] += vs_num
        end
    end
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        buffer[1] = MPI.Reduce(buffer[1],
            MPI.MAX,
            MPI.COMM_WORLD
        )
        buffer[2] = MPI.Reduce(buffer[1],
            +,
            MPI.COMM_WORLD
        )
    else
        MPI.Reduce(buffer[1],
            MPI.MAX,
            MPI.COMM_WORLD
        )
        MPI.Reduce(buffer[2],
            +,
            MPI.COMM_WORLD
        )
    end
    return buffer
end