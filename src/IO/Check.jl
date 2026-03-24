"""
$(SIGNATURES)
Perform a check procedure. Simulation status will be output every `ST_CHECK_INTERVAL` steps.
"""
function check!(i,ps4est,amr)
    if amr.global_data.status.residual.step%amr.global_data.config.solver.ST_CHECK_INTERVAL==0
        max_vs_num,total_phase_num = check_vs_num(amr)
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            i+=1
            println("Iteration: $i")
            sim_time = amr.global_data.status.sim_time
            println("Simulation time: $sim_time")
            res = maximum(amr.global_data.status.residual.residual)
            println("Residual: $res")
            ref_vs_num = amr.global_data.status.max_vs_num
            println("MPI buffer size: $ref_vs_num")
            println("Maximum number of velocity grids: $max_vs_num")
            pp = PointerWrapper(ps4est)
            global_num_quadrants = pp.global_num_quadrants[]
            println("Total number of physical grids: $global_num_quadrants")
            println("Total number of phase grids: $total_phase_num")
        end
        check_for_save!(ps4est,amr)
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
function check_for_animsave!(ps4est::Ptr{p4est_t},amr;path="./animation")
    sim_time = amr.global_data.status.sim_time
    output = amr.global_data.config.output
    if sim_time==0.
        vertices,cells,point_solutions,solutions = pvtu_data(ps4est,amr,output.vtk_celltype)
        ranks = ones(Int,size(solutions,1))*MPI.Comm_rank(MPI.COMM_WORLD)
        step = output.anim_index
        if MPI.Comm_rank(MPI.COMM_WORLD)==0
            paraview_collection(path*"/full_simulation";append=step!=0) do pvd
                pvtk_grid(path*"/step$step",vertices,cells;part = MPI.Comm_rank(MPI.COMM_WORLD)+1,nparts = MPI.Comm_size(MPI.COMM_WORLD)) do pvtk
                    pvtk["rho"] = @views solutions[:,1]
                    pvtk["velocity"] = @views (solutions[:,2],solutions[:,3])
                    pvtk["T"] = @views solutions[:,4]
                    pvtk["qf"] = (solutions[:,5],solutions[:,6])
                    pvtk["mpi_rank"] = ranks
                    pvtk["rho",VTKPointData()] = @views point_solutions[:,1]
                    pvtk["velocity",VTKPointData()] = @views (point_solutions[:,2],point_solutions[:,3])
                    pvtk["T",VTKPointData()] = @views point_solutions[:,4]
                    pvtk["qf",VTKPointData()] = (point_solutions[:,5],point_solutions[:,6])
                    pvd[sim_time] = pvtk
                    close(pvd)
                end
            end
        else
            pvtk_grid(path*"/step$step",vertices,cells;part = MPI.Comm_rank(MPI.COMM_WORLD)+1,nparts = MPI.Comm_size(MPI.COMM_WORLD)) do pvtk
                pvtk["rho"] = @views solutions[:,1]
                pvtk["velocity"] = @views (solutions[:,2],solutions[:,3])
                pvtk["T"] = @views solutions[:,4]
                pvtk["qf"] = (solutions[:,5],solutions[:,6])
                pvtk["mpi_rank"] = ranks
                pvtk["rho",VTKPointData()] = @views point_solutions[:,1]
                pvtk["velocity",VTKPointData()] = @views (point_solutions[:,2],point_solutions[:,3])
                pvtk["T",VTKPointData()] = @views point_solutions[:,4]
                pvtk["qf",VTKPointData()] = (point_solutions[:,5],point_solutions[:,6])
            end
        end
        if output.vs_output_criterion!=null_udf
            trees = amr.field.trees.data
            for tree in trees
                for ps_data in tree
                    (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
                    id,flag = output.vs_output_criterion(;ps_data,amr)
                    if flag
                        vertices,cells,point_solutions,solutions = vtk_data(ps_data.vs_data,amr,output.vs_vtk_celltype)
                        paraview_collection(path*"/id$(id)vs";append=step!=0) do pvd
                            vtk_grid(path*"/id$(id)step$(step)vs",vertices,cells) do vtk
                                vtk["df"] = @views Tuple([v for v in eachcol(solutions)])
                                vtk["df",VTKPointData()] = @views Tuple([v for v in eachcol(point_solutions)])
                                pvd[sim_time] = vtk
                                close(pvd)
                            end
                        end
                    end
                end
            end
        end
    elseif div(sim_time,output.anim_dt)==output.anim_index+1
        vertices,cells,point_solutions,solutions = pvtu_data(ps4est,amr,output.vtk_celltype)
        ranks = ones(Int,size(solutions,1))*MPI.Comm_rank(MPI.COMM_WORLD)
        step = output.anim_index+1
        if MPI.Comm_rank(MPI.COMM_WORLD)==0
            paraview_collection(path*"/full_simulation";append=step!=1) do pvd
                pvtk_grid(path*"/step$step",vertices,cells;part = MPI.Comm_rank(MPI.COMM_WORLD)+1,nparts = MPI.Comm_size(MPI.COMM_WORLD)) do pvtk
                    pvtk["rho"] = @views solutions[:,1]
                    pvtk["velocity"] = @views (solutions[:,2],solutions[:,3])
                    pvtk["T"] = @views solutions[:,4]
                    pvtk["qf"] = (solutions[:,5],solutions[:,6])
                    pvtk["mpi_rank"] = ranks
                    pvtk["rho",VTKPointData()] = @views point_solutions[:,1]
                    pvtk["velocity",VTKPointData()] = @views (point_solutions[:,2],point_solutions[:,3])
                    pvtk["T",VTKPointData()] = @views point_solutions[:,4]
                    pvtk["qf",VTKPointData()] = (point_solutions[:,5],point_solutions[:,6])
                    pvd[sim_time] = pvtk
                    close(pvd)
                end
            end
        else
            pvtk_grid(path*"/step$step",vertices,cells;part = MPI.Comm_rank(MPI.COMM_WORLD)+1,nparts = MPI.Comm_size(MPI.COMM_WORLD)) do pvtk
                pvtk["rho"] = @views solutions[:,1]
                pvtk["velocity"] = @views (solutions[:,2],solutions[:,3])
                pvtk["T"] = @views solutions[:,4]
                pvtk["qf"] = (solutions[:,5],solutions[:,6])
                pvtk["mpi_rank"] = ranks
                pvtk["rho",VTKPointData()] = @views point_solutions[:,1]
                pvtk["velocity",VTKPointData()] = @views (point_solutions[:,2],point_solutions[:,3])
                pvtk["T",VTKPointData()] = @views point_solutions[:,4]
                pvtk["qf",VTKPointData()] = (point_solutions[:,5],point_solutions[:,6])
            end
        end
        if output.vs_output_criterion!=null_udf
            trees = amr.field.trees.data
            for tree in trees
                for ps_data in tree
                    (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
                    id,flag = output.vs_output_criterion(;ps_data,amr)
                    if flag
                        vertices,cells,point_solutions,solutions = vtk_data(ps_data.vs_data,amr,output.vs_vtk_celltype)
                        paraview_collection(path*"/id$(id)vs";append=step!=1) do pvd
                            vtk_grid(path*"/id$(id)step$(step)vs",vertices,cells) do vtk
                                vtk["df"] = @views Tuple([v for v in eachcol(solutions)])
                                vtk["df",VTKPointData()] = @views Tuple([v for v in eachcol(point_solutions)])
                                pvd[sim_time] = vtk
                                close(pvd)
                            end
                        end
                    end
                end
            end
        end
        output.anim_index += 1
    end
end

function check_vs_num(amr::KitAMR_Data)
    trees = amr.field.trees.data
    buffer = zeros(Int,2) # [max,total]
    for tree in trees
        for ps_data in tree
            isa(ps_data,InsideSolidData)&&continue
            ps_data.bound_enc<0&&continue
            vs_num = ps_data.vs_data.vs_num
            buffer[1] = max(buffer[1],vs_num)
            buffer[2] += vs_num
        end
    end
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        buffer[1] = MPI.Reduce(buffer[1],(x,y)->(
            max(x,y)
        ),
        MPI.COMM_WORLD
        )
        buffer[2] = MPI.Reduce(buffer[1],(x,y)->(
            x+y
        ),
        MPI.COMM_WORLD
        )
    else
        MPI.Reduce(buffer[1],(x,y)->(
            max(x,y)
        ),
        MPI.COMM_WORLD
        )
        MPI.Reduce(buffer[2],(x,y)->(
            x+y
        ),
        MPI.COMM_WORLD
        )
    end
    return buffer
end