function read_config(filename::T) where {T<:AbstractString}
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @info "reading config from $filename"
        println("--------------------------------------------------------------")
    end
    f = open(filename,"r";lock = false)
    config_dict = Dict{Symbol,Any}()
    for line in eachline(f)
        if length(line) == 0 || line[1] == '#'
            continue
        end
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            println(line)
        end
        var, val = split(line, "=")
        stripped = Symbol(strip(var))
        val = eval(Meta.parse(split(val, "#")[1]))
        config_dict[stripped] = val
    end
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        println("--------------------------------------------------------------")
        println("")
    end
    close(f)
    return config_dict
end

function string_to_Cstring(s::String)
    t = codeunits(s)
    p = Ptr{Cstring}(malloc(sizeof(UInt8)*(length(t)+1)))
    ap = Base.unsafe_wrap(Vector{UInt8},Ptr{UInt8}(p),length(t)+1)
    ap[1:end-1] .= t
    ap[end] = 0
    return p
end

#=
example:
function fieldvalues_fn(ps_data)
    return [1/ps_data.prim[end]]
end

write_VTK(ps4est,"testT",["T"],fieldvalues_fn)
=#
function write_VTK(ps4est::Ptr{p4est_t},filename::String,fieldnames::Vector{String},fieldvalues_fn)
    p4est_geo = p4est_geometry_new_connectivity(pointer(PointerWrapper(ps4est).connectivity))
    cont = p4est_vtk_context_new(ps4est,filename)
    p4est_vtk_context_set_geom(cont,p4est_geo)
    p4est_vtk_context_set_continuous(cont,1)
    cont = p4est_vtk_write_header(cont)
    fp = PointerWrapper(ps4est)
    num_quads = fp.local_num_quadrants[]
    ps = Vector{Ptr{Nothing}}(undef,length(fieldnames))
    pas = Vector{Vector{Cdouble}}(undef,length(fieldnames))
    pscs = Vector{Ptr{sc_array_t}}(undef,length(fieldnames))
    for i in eachindex(fieldnames)
        ps[i] = sc_malloc(-1,num_quads*8)
        pas[i] = unsafe_wrap(Vector{Cdouble},Ptr{Cdouble}(ps[i]),num_quads)
        pscs[i] = sc_array_new_data(ps[i],8,num_quads)
    end
    function init_cell_data_kernel(ip,data,dp)
        pas = unsafe_pointer_to_objref(data)
        qid = local_quadid(ip)+1
        ps_data = Base.unsafe_pointer_to_objref(pointer(dp.ps_data))
        fieldvalues = fieldvalues_fn(ps_data)
        for i in eachindex(fieldvalues)
            pas[i][qid] = fieldvalues[i]
        end
    end
    function init_cell_data(info,data)
        AMR_volume_iterate(info, data, P4est_PS_Data, init_cell_data_kernel)
    end
    ppas = pointer_from_objref(pas)
    cfn = @cfunction($init_cell_data, Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})) # explicit cfn to avoid closure error on Windows and MacOS platforms.
    GC.@preserve ps4est ppas pas p4est_iterate(
        ps4est,
        C_NULL,
        ppas,
        cfn,
        C_NULL,
        C_NULL,
    )
    pairs = Vector{Union{Ptr{Cstring},Ptr{sc_array_t}}}(undef,2*length(fieldnames))
    p_names = Vector{Ptr{Cstring}}(undef,length(fieldnames))
    for i in eachindex(fieldnames)
        pairs[2*i-1] = p_names[i] = string_to_Cstring(fieldnames[i])
        pairs[2*i] = pscs[i]
    end
    cont = p4est_vtk_write_cell_dataf(cont,1,1,1,0,length(fieldnames),0,pairs...,cont)
    p4est_vtk_write_footer(cont)
    for i in eachindex(fieldnames)
        sc_free(-1,Ptr{Nothing}(p_names[i]))
        sc_free(-1,ps[i])
    end
end
function write_VTK(ps4est::Ptr{p8est_t},filename::String,fieldnames::Vector{String},fieldvalues_fn)
    p4est_geo = p8est_geometry_new_connectivity(pointer(PointerWrapper(ps4est).connectivity))
    cont = p8est_vtk_context_new(ps4est,filename)
    p8est_vtk_context_set_geom(cont,p4est_geo)
    p8est_vtk_context_set_continuous(cont,1)
    cont = p8est_vtk_write_header(cont)
    fp = PointerWrapper(ps4est)
    num_quads = fp.local_num_quadrants[]
    ps = Vector{Ptr{Nothing}}(undef,length(fieldnames))
    pas = Vector{Vector{Cdouble}}(undef,length(fieldnames))
    pscs = Vector{Ptr{sc_array_t}}(undef,length(fieldnames))
    for i in eachindex(fieldnames)
        ps[i] = sc_malloc(-1,num_quads*8)
        pas[i] = unsafe_wrap(Vector{Cdouble},Ptr{Cdouble}(ps[i]),num_quads)
        pscs[i] = sc_array_new_data(ps[i],8,num_quads)
    end
    function init_cell_data_kernel(ip,data,dp)
        pas = unsafe_pointer_to_objref(data)
        qid = local_quadid(ip)+1
        ps_data = Base.unsafe_pointer_to_objref(pointer(dp.ps_data))
        fieldvalues = fieldvalues_fn(ps_data)
        for i in eachindex(fieldvalues)
            pas[i][qid] = fieldvalues[i]
        end
    end
    function init_cell_data(info,data)
        AMR_volume_iterate(info, data, P4est_PS_Data, init_cell_data_kernel)
    end
    ppas = pointer_from_objref(pas)
    cfn = @cfunction($init_cell_data, Cvoid, (Ptr{p8est_iter_volume_info}, Ptr{Nothing})) # explicit cfn to avoid closure error on Windows and MacOS platforms.
    GC.@preserve ps4est ppas pas p8est_iterate(
        ps4est,
        C_NULL,
        ppas,
        cfn,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    pairs = Vector{Union{Ptr{Cstring},Ptr{sc_array_t}}}(undef,2*length(fieldnames))
    p_names = Vector{Ptr{Cstring}}(undef,length(fieldnames))
    for i in eachindex(fieldnames)
        pairs[2*i-1] = p_names[i] = string_to_Cstring(fieldnames[i])
        pairs[2*i] = pscs[i]
    end
    cont = p8est_vtk_write_cell_dataf(cont,1,1,1,0,length(fieldnames),0,pairs...,cont)
    p8est_vtk_write_footer(cont)
    for i in eachindex(fieldnames)
        sc_free(-1,Ptr{Nothing}(p_names[i]))
        sc_free(-1,ps[i])
    end
end

#=
example:
function fieldvalues_fn(vs_data::VS_Data{3})
    return [vs_data.df[:,1],vs_data.level]
end
=#
function write_vs_VTK(vs_data::AbstractVsData{2,2},amr::AMR{2,2},filename::String,fieldnames::Vector{String},fieldvalues_fn)
    global_data = amr.global_data
    xmin,xmax,ymin,ymax = global_data.config.quadrature
    Nx,Ny = global_data.config.vs_trees_num
    AMR_VS_MAXLEVEL = global_data.config.solver.AMR_VS_MAXLEVEL
    dx = (xmax - xmin) / Nx/2^AMR_VS_MAXLEVEL
    dy = (ymax - ymin) / Ny/2^AMR_VS_MAXLEVEL
    D = [dx,dy]
    vertices = Matrix{Float64}(undef,2,8*length(vs_data.level))
    midpoint = vs_data.midpoint
    level = vs_data.level
    dlevel = -(vs_data.level .- AMR_VS_MAXLEVEL)
    cells = Vector{MeshCell}(undef,length(level))
    for i in eachindex(level)
        for j in 1:4
            @. vertices[:,(i-1)*4+j] = midpoint[i,:]+RMT[2][j]*2^dlevel[i]*D/2
        end
        cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL,(1:4).+4*(i-1))
    end
    vtk_grid(filename,vertices,cells;append=false) do vtk
        cell_datas = fieldvalues_fn(vs_data)
        for i in eachindex(fieldnames)
            vtk[fieldnames[i],VTKCellData()] = cell_datas[i]
        end
    end
end
function write_vs_VTK(vs_data::AbstractVsData{3,1},amr::AMR{3,1},filename::String,fieldnames::Vector{String},fieldvalues_fn)
    global_data = amr.global_data
    xmin,xmax,ymin,ymax,zmin,zmax = global_data.config.quadrature
    Nx,Ny,Nz = global_data.config.vs_trees_num
    AMR_VS_MAXLEVEL = global_data.config.solver.AMR_VS_MAXLEVEL
    dx = (xmax - xmin) / Nx/2^AMR_VS_MAXLEVEL
    dy = (ymax - ymin) / Ny/2^AMR_VS_MAXLEVEL
    dz = (zmax - zmin) / Nz/2^AMR_VS_MAXLEVEL
    D = [dx,dy,dz]
    vertices = Matrix{Float64}(undef,3,8*length(vs_data.level))
    midpoint = vs_data.midpoint
    level = vs_data.level
    dlevel = -(vs_data.level .- AMR_VS_MAXLEVEL)
    cells = Vector{MeshCell}(undef,length(level))
    for i in eachindex(level)
        for j in 1:8
            @. vertices[:,(i-1)*8+j] = midpoint[i,:]+RMT[3][j]*2^dlevel[i]*D/2
        end
        cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL,(1:8).+8*(i-1))
    end
    vtk_grid(filename,vertices,cells;append=false) do vtk
        cell_datas = fieldvalues_fn(vs_data)
        for i in eachindex(fieldnames)
            vtk[fieldnames[i],VTKCellData()] = cell_datas[i]
        end
    end
end
function neighbor_num(ps_data::PS_Data,::P_pxest_t,::AMR,::Integer)
    return abs.(ps_data.neighbor.state)
end
function neighbor_num(::InsideSolidData,ps4est::P_pxest_t,amr::AMR{DIM},quadid::Integer) where {DIM}
    neighbor_num = Vector{Int}(undef,2*DIM)
    global_data = amr.global_data
    ghost = global_data.forest.ghost
    mesh = global_data.forest.mesh
    for dir = 1:2*DIM
        neighbor_quads = sc_array_new(sizeof(P_pxest_quadrant_t))
        neighbor_encs = sc_array_new(sizeof(Cint))
        neighbor_qid = sc_array_new(sizeof(Cint))
        GC.@preserve dir ps4est ghost mesh quadid neighbor_encs neighbor_qid neighbor_quads p4est_mesh_get_neighbors(
            ps4est,
            ghost,
            mesh,
            quadid,
            dir-1,
            neighbor_quads,
            neighbor_encs,
            neighbor_qid,
        )
        @inbounds neighbor_num[dir] = PointerWrapper(neighbor_encs).elem_count[]
        sc_array_destroy(neighbor_quads)
        sc_array_destroy(neighbor_encs)
        sc_array_destroy(neighbor_qid)
    end
    return neighbor_num
end
function save_result(ps4est::P_pxest_t,amr::AMR{DIM,NDF}) where{DIM,NDF}
    fp = PointerWrapper(ps4est)
    ps_solution = Vector{PS_Solution}(undef,fp.local_num_quadrants[])
    neighbor_nums = Vector{Vector{Int}}(undef,fp.local_num_quadrants[])
    trees = amr.field.trees.data
    config = amr.global_data.config
    index = 1
    for i in eachindex(trees)
        for j in eachindex(trees[i])
            ps_data = trees[i][j]
            ps_solution[index] = PS_Solution(ps_data)
            neighbor_nums[index] = neighbor_num(ps_data,ps4est,amr,index-1)
            index+=1
        end
    end
    solution = Solution(ps_solution)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    result = Result(solution,MeshInfo(neighbor_nums))
    MPI.Barrier(MPI.COMM_WORLD)
    dir_path = "./result"*Dates.format(now(), "yyyy-mm-dd_HH-MM")*"/"
    !isdir(dir_path) && mkpath(dir_path)
    p4est_save_ext("p",ps4est,Cint(0),Cint(0))
    if rank==0
        size = MPI.Comm_size(MPI.COMM_WORLD)
        solverset = SolverSet(ConfigureForSave(config),size)
        save_object(dir_path * "solverset.jld2", solverset)
    end
    save_object(dir_path * "result_"*string(rank)*".jld2", result)
    save_boundary_result(dir_path,amr)
end
# function save_boundary_result(dir_path::String,amr::AMR{DIM,NDF}) where{DIM,NDF}
#     boundary = amr.field.boundary
#     ibs = amr.global_data.config.IB
#     boundary_results = Vector{Boundary_Solution}(undef,length(ibs))
#     for i in eachindex(ibs)
#         ib = ibs[i]
#         solid_cells = boundary.solid_cells[i]
#         ib_nodes = boundary.IB_cells[i].IB_nodes
#         midpoints = Vector{Vector{Float64}}(undef,length(solid_cells.ps_datas))
#         ps_solutions = Vector{Boundary_PS_Solution}(undef,length(solid_cells.ps_datas))
#         for j in eachindex(solid_cells.ps_datas)
#             midpoints[j],ps_solutions[j] = save_boundary_result(ib,solid_cells.ps_datas[j],ib_nodes[j],boundary.IB_cells[i].templates[j],amr.global_data)
#         end
#         boundary_results[i] = Boundary_Solution(midpoints,ps_solutions)
#     end
#     rank = MPI.Comm_rank(MPI.COMM_WORLD)
#     save_object(dir_path*"boundary_result_"*string(rank)*".jld2",boundary_results)
# end
function save_boundary_result(dir_path::String,amr::AMR{DIM,NDF}) where{DIM,NDF}
    ibs = amr.global_data.config.IB
    boundary_results = [Boundary_Solution(Vector{Float64}[],Boundary_PS_Solution[])]
    for tree in amr.field.trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<=0)&&continue
            ib = ibs[ps_data.bound_enc]
            save_boundary_result!(ib,ps_data,boundary_results,amr)
        end
    end
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    save_object(dir_path*"boundary_result_"*string(rank)*".jld2",boundary_results)
end
function boundary_write_csv(csvname,results,config::ConfigureForSave{2})
    for i in eachindex(config.IB)
        df = DataFrame()
        df.x=[x[1] for x in results[i].midpoints];df.y=[x[2] for x in results[i].midpoints]
        df.rho=[x.prim[1] for x in results[i].ps_solutions]
        df.u=[x.prim[2] for x in results[i].ps_solutions]
        df.v=[x.prim[3] for x in results[i].ps_solutions]
        df.T=[1/x.prim[4] for x in results[i].ps_solutions]
        df.qfx=[x.qf[1] for x in results[i].ps_solutions];df.qfy=[x.qf[2] for x in results[i].ps_solutions]
        df.p11 = [x.p[1] for x in results[i].ps_solutions];df.p12 = [x.p[2] for x in results[i].ps_solutions]
        df.p22 = [x.p[3] for x in results[i].ps_solutions]
        CSV.write(csvname*"_"*string(i)*".csv",df)
    end
end
function boundary_result2csv(dirname::String,csvname::String)
    path = "./"*dirname
    solverset = load_object(path*"/solverset.jld2")
    results = nothing
    for i in 1:solverset.mpi_size
        if i==1
            results = load_object(path*"/boundary_result_"*string(i-1)*".jld2")
        else
            rs = load_object(path*"/boundary_result_"*string(i-1)*".jld2")
            for j in eachindex(rs)
                r = rs[j]
                result = results[j]
                append!(result.ps_solutions,r.ps_solutions)
                append!(result.midpoints,r.midpoints)
            end
        end
    end
    boundary_write_csv(csvname,results,solverset.config)
end
function result2vtk(dirname::String,vtkname::String)
    !MPI.Initialized() && MPI.Init()
    path = "./"*dirname
    solverset = load_object(path*"/solverset.jld2")
    if typeof(solverset.config).parameters[1]==2
        DIM=2
        cnn = Ptr{Ptr{p4est_connectivity_t}}(Libc.malloc(sizeof(Ptr{Ptr{p4est_connectivity_t}})))
        ps4est = p4est_load_ext("p",MPI.COMM_WORLD,Cint(0),Cint(0),Cint(1),Cint(0),C_NULL,cnn)
        result = nothing
        ranks = nothing
        for i = 1:solverset.mpi_size
            if i==1
                result = load_object(path*"/result_"*string(i-1)*".jld2")
                ranks = zeros(Int,length(result.solution.ps_solutions))
            else
                r = load_object(path*"/result_"*string(i-1)*".jld2")
				append!(result.solution.ps_solutions,r.solution.ps_solutions)
				append!(result.solution.vs_solutions,r.solution.vs_solutions)
				append!(result.mesh_info.neighbor_nums,r.mesh_info.neighbor_nums)
                append!(ranks,ones(Int,length(r.solution.ps_solutions))*(i-1))
            end
        end
        vtk_cnn = Vector{Vector{Int}}(undef,length(result.solution.ps_solutions))
        neighbor_nums = result.mesh_info.neighbor_nums
        for i in eachindex(neighbor_nums)
            for j in eachindex(neighbor_nums[i])
                neighbor_nums[i][j]==0&&(neighbor_nums[i][j]=1)
            end
        end
        for i in eachindex(vtk_cnn)
            addi_points_num = 0
            for j in eachindex(neighbor_nums[i])
                addi_points_num+= neighbor_nums[i][j]==2^(DIM-1) ? 1 : 0
            end
            vtk_cnn[i] = Vector{Int}(undef,2^DIM+addi_points_num)
        end
        points = Vector{Float64}[]
        data = Vector{Any}()
        for el in (points,vtk_cnn,neighbor_nums)
            push!(data,el)
        end
        p_data = pointer_from_objref(data)
        GC.@preserve data AMR_corner_iterate(ps4est;user_data = p_data) do ip,data
            points,vtk_cnn,neighbor_nums = unsafe_pointer_to_objref(data)
            DIM=isa(ip,PointerWrapper{p4est_iter_corner_info_t}) ? 2 : 3
            for i in 1:ip.sides.elem_count[]
                side = iPointerWrapper(ip.sides,p4est_iter_corner_side_t,i-1)
                cornerid = side.corner[]+1 # z-order
                if i==1
                    ds,midpoint = quad_to_cell(ip.p4est,side.treeid[],side.quad)
                    point = @. midpoint+0.5*ds*RMT[DIM][cornerid]
                    push!(points,point)
                end
                id = length(points)
                quadid = local_quadid(ip,side)
                neighbor_num = neighbor_nums[quadid+1]
                cornerid==1&&(vtk_cnn[quadid+1][1]=id)
                cornerid==2&&(vtk_cnn[quadid+1][neighbor_num[3]+1]=id)
                cornerid==3&&(vtk_cnn[quadid+1][sum(neighbor_num[2:4])+1]=id)
                cornerid==4&&(vtk_cnn[quadid+1][sum(neighbor_num[2:3])+1]=id)
            end
            return nothing
        end
        GC.@preserve data AMR_face_iterate(ps4est;user_data = p_data) do ip,data
            points,vtk_cnn,neighbor_nums = unsafe_pointer_to_objref(data)
            DIM=isa(ip,PointerWrapper{p4est_iter_face_info_t}) ? 2 : 3
            ip.sides.elem_count[]==1&&return nothing
            side1 = iPointerWrapper(ip.sides,p4est_iter_face_side_t,0)
            side2 = iPointerWrapper(ip.sides,p4est_iter_face_side_t,1)
            (side1.is_hanging[]==0&&side2.is_hanging[]==0)&&return nothing
            for side in (side1,side2)
                    if side.is_hanging[]==0
                        ds,midpoint = quad_to_cell(ip.p4est,side.treeid[],side.is.full.quad)
                        faceid = side.face[]+1
                        point = @. midpoint+0.5*ds*NMT[DIM][faceid]
                        push!(points,point)
                        id = length(points)
                        quadid = local_quadid(ip.p4est,side.treeid[],side.is.full.quadid[])
                        neighbor_num = neighbor_nums[quadid+1]
                        faceid==1&&(@inbounds vtk_cnn[quadid+1][end]=id)
                        faceid==2&&(@inbounds vtk_cnn[quadid+1][neighbor_num[3]+2]=id)
                        faceid==3&&(@inbounds vtk_cnn[quadid+1][2]=id)
                        faceid==4&&(@inbounds vtk_cnn[quadid+1][sum(neighbor_num[2:3])+2]=id)
                    end
            end
            for side in (side1,side2)
                    if side.is_hanging[]==1
                        faceid = side.face[]+1
                        id = length(points)
                        quadids = side.is.hanging.quadid[]
                        quadid1 = local_quadid(ip.p4est,side.treeid[],quadids[1])
                        quadid2 = local_quadid(ip.p4est,side.treeid[],quadids[2])
                        neighbor_num1 = neighbor_nums[quadid1+1]
                        neighbor_num2 = neighbor_nums[quadid2+1]
                        faceid==1&&(@inbounds vtk_cnn[quadid1+1][end]=id;vtk_cnn[quadid2+1][1]=id)
                        faceid==2&&(@inbounds vtk_cnn[quadid1+1][sum(neighbor_num1[2:3])+1]=id;vtk_cnn[quadid2+1][neighbor_num2[3]+1]=id)
                        faceid==3&&(@inbounds vtk_cnn[quadid1+1][2]=id;vtk_cnn[quadid2+1][1]=id)
                        faceid==4&&(@inbounds vtk_cnn[quadid1+1][sum(neighbor_num1[2:3])+1]=id;vtk_cnn[quadid2+1][end]=id)
                    end
            end
            return nothing
        end
        cells = [MeshCell(PolyData.Polys(),cnn) for cnn in vtk_cnn]
        vertices = Matrix{Float64}(undef,2,length(points))
        for i in eachindex(points)
            @inbounds vertices[:,i] .= points[i]
        end
        vtk_grid(vtkname,vertices,cells) do vtk
            vtk["rho"] = [ps_solution.prim[1] for ps_solution in result.solution.ps_solutions]
            # vtk["u"] = [ps_solution.prim[2] for ps_solution in result.solution.ps_solutions]
            # vtk["v"] = [ps_solution.prim[3] for ps_solution in result.solution.ps_solutions]
            vtk["velocity"] = ([ps_solution.prim[2] for ps_solution in result.solution.ps_solutions],
                [ps_solution.prim[3] for ps_solution in result.solution.ps_solutions])
            vtk["T"] = [1/ps_solution.prim[end] for ps_solution in result.solution.ps_solutions]
            # vtk["qfx"] = [ps_solution.qf[1] for ps_solution in result.solution.ps_solutions]
            # vtk["qfy"] = [ps_solution.qf[2] for ps_solution in result.solution.ps_solutions]
            vtk["qf"] = ([ps_solution.qf[1] for ps_solution in result.solution.ps_solutions],
                [ps_solution.qf[2] for ps_solution in result.solution.ps_solutions])
            vtk["mpi_rank"] = ranks
        end
		fp = PointerWrapper(ps4est)
        p4est_connectivity_destroy(pointer(fp.connectivity))
        p4est_destroy(ps4est)
    else
        @error "Only support 2D now"
    end
end