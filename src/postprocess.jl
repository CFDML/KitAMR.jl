function collect_results(ps4est, DVM_data)
    pp = PointerWrapper(ps4est)
    trees = DVM_data.trees
    solutions = Array{Solution}(undef, pp.local_num_quadrants[])
    index = 1
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            solutions[index] = Solution(
                SVector{DIM,Cdouble}(ps_data.midpoint),
                SVector{DIM + 2,Cdouble}(ps_data.w),
                SVector{DIM + 2,Cdouble}(ps_data.prim),
                SVector{DIM,Cdouble}(ps_data.qf),
            )
            index += 1
        end
    end
    return solutions
end
function communicate_solution(solutions::Vector{T}, num_quadrants) where {T}
    if MPI.Comm_rank(MPI.COMM_WORLD) != 0
        buffer = solutions
        sreq =
            MPI.Isend(buffer, MPI.COMM_WORLD; dest = 0, tag = MPI.Comm_rank(MPI.COMM_WORLD))
        stats = MPI.Wait(sreq)
        return nothing
    else
        for i = 2:MPI.Comm_size(MPI.COMM_WORLD)
            buffer = Array{T}(undef, num_quadrants[i])
            rreq = MPI.Irecv!(buffer, MPI.COMM_WORLD; source = i - 1, tag = i - 1)
            state = MPI.Wait(rreq)
            append!(solutions, buffer)
        end
        return solutions
    end
end
function collect_num_quadrants(ps4est::Ptr{p4est_t})
    pp = PointerWrapper(ps4est)
    if MPI.Comm_rank(MPI.COMM_WORLD) != 0
        buffer = Array{Cint}(undef, 1)
        buffer[1] = pp.local_num_quadrants[]
        sreq =
            MPI.Isend(buffer, MPI.COMM_WORLD; dest = 0, tag = MPI.Comm_rank(MPI.COMM_WORLD))
        stats = MPI.Wait(sreq)
        return nothing
    else
        num_quadrants = zeros(Cint, MPI.Comm_size(MPI.COMM_WORLD))
        buffer = Array{Cint}(undef, 1)
        for i = 2:MPI.Comm_size(MPI.COMM_WORLD)
            rreq = MPI.Irecv!(buffer, MPI.COMM_WORLD; source = i - 1, tag = i - 1)
            state = MPI.Wait(rreq)
            num_quadrants[i] = buffer[1]
        end
        return num_quadrants
    end
end
function collect_solution(ps4est, DVM_data)
    num_quadrants = collect_num_quadrants(ps4est)
    solutions = collect_results(ps4est, DVM_data)
    solutions = communicate_solution(solutions, num_quadrants)
    return solutions
end
function reshape_solutions(
    solutions::Vector,
    global_data::Global_Data,
    field::Symbol,
    index::Int,
)
    Nx, Ny = global_data.trees_num
    xmin, xmax, ymin, ymax = global_data.geometry
    midpoints = zeros(length(solutions), 2)
    z = zeros(length(solutions))
    for i in eachindex(solutions)
        midpoints[i, :] .= solutions[i].midpoint
        z[i] = getfield(solutions[i], field)[index]
    end
    itp = scipy.interpolate.CloughTocher2DInterpolator(midpoints, z)
    dx = (xmax - xmin) / Nx / 2^DVM_PS_MAXLEVEL
    dy = (ymax - ymin) / Ny / 2^DVM_PS_MAXLEVEL
    X = collect(xmin+EPS+dx/2:dx:xmax-EPS-dx/2)
    Y = collect(ymin+EPS+dy/2:dy:ymax-EPS-dy/2)
    pyZ = np.zeros((length(X), length(Y)))
    for i in eachindex(X)
        for j in eachindex(Y)
            pyZ[i-1, j-1] = itp(X[i], Y[j])[]
        end
    end
    Z = pyconvert(Matrix{Float64}, pyZ)
    return X, Y, Z
end

function collect_mesh(DVM_data::DVM_Data)
    mesh_points = collect_mesh_points(DVM_data)
    num_faces = collect_num_faces(mesh_points)
    meshes = communicate_mesh(mesh_points, num_faces)
    return meshes
end
function mesh_plot(x::AV, y::AV, ax::Axis)
    for i = 1:Int(length(x) / 2)
        # x1 = x[2*i-1];x2 = x[2*i]
        # y1 = y[2*i-1];y2 = y[2*i]
        lines!(ax, @view(x[2*i-1:2*i]), @view(y[2*i-1:2*i]), color = :red, linewidth = 0.5)
    end
end
function communicate_mesh(mesh_points::T, num_faces) where {T}
    if MPI.Comm_rank(MPI.COMM_WORLD) != 0
        buffer = mesh_points
        sreq =
            MPI.Isend(buffer, MPI.COMM_WORLD; dest = 0, tag = MPI.Comm_rank(MPI.COMM_WORLD))
        stats = MPI.Wait(sreq)
        return nothing
    else
        meshes = Array{Mesh_Points}(undef, 1)
        meshes[1] = mesh_points
        for i = 2:MPI.Comm_size(MPI.COMM_WORLD)
            buffer = Array{Mesh_Points{SVector{num_faces[i],Float64}}}(undef, 1)
            rreq = MPI.Irecv!(buffer, MPI.COMM_WORLD; source = i - 1, tag = i - 1)
            state = MPI.Wait(rreq)
            append!(meshes, buffer)
        end
        return meshes
    end
end
function reshape_mesh(meshes::Vector{Mesh_Points})
    x = Vector{Float64}(undef, 0)
    y = Vector{Float64}(undef, 0)
    for i in eachindex(meshes)
        append!(x, collect(meshes[i].x))
        append!(y, collect(meshes[i].y))
    end
    return x, y
end
function collect_num_faces(mesh_points::Mesh_Points{SVector{N,T}}) where {N,T}
    if MPI.Comm_rank(MPI.COMM_WORLD) != 0
        buffer = Array{Int}(undef, 1)
        buffer[1] = N
        sreq =
            MPI.Isend(buffer, MPI.COMM_WORLD; dest = 0, tag = MPI.Comm_rank(MPI.COMM_WORLD))
        stats = MPI.Wait(sreq)
        return nothing
    else
        num_faces = zeros(Int, MPI.Comm_size(MPI.COMM_WORLD))
        buffer = Array{Int}(undef, 1)
        for i = 2:MPI.Comm_size(MPI.COMM_WORLD)
            rreq = MPI.Irecv!(buffer, MPI.COMM_WORLD; source = i - 1, tag = i - 1)
            state = MPI.Wait(rreq)
            num_faces[i] = buffer[1]
        end
        return num_faces
    end
end
function collect_mesh_points_kernal!(::Val{1}, ps_data::AbstractPsData, x::AV, y::AV)
    x1 = ps_data.midpoint[1] - 0.5 * ps_data.ds[1]
    x2 = ps_data.midpoint[1] - 0.5 * ps_data.ds[1]
    y1 = ps_data.midpoint[2] - 0.5 * ps_data.ds[2]
    y2 = ps_data.midpoint[2] + 0.5 * ps_data.ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_mesh_points_kernal!(::Val{2}, ps_data::AbstractPsData, x::AV, y::AV)
    x1 = ps_data.midpoint[1] + 0.5 * ps_data.ds[1]
    x2 = ps_data.midpoint[1] + 0.5 * ps_data.ds[1]
    y1 = ps_data.midpoint[2] - 0.5 * ps_data.ds[2]
    y2 = ps_data.midpoint[2] + 0.5 * ps_data.ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_mesh_points_kernal!(::Val{3}, ps_data::AbstractPsData, x::AV, y::AV)
    x1 = ps_data.midpoint[1] - 0.5 * ps_data.ds[1]
    x2 = ps_data.midpoint[1] + 0.5 * ps_data.ds[1]
    y1 = ps_data.midpoint[2] - 0.5 * ps_data.ds[2]
    y2 = ps_data.midpoint[2] - 0.5 * ps_data.ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_mesh_points_kernal!(::Val{4}, ps_data::AbstractPsData, x::AV, y::AV)
    x1 = ps_data.midpoint[1] - 0.5 * ps_data.ds[1]
    x2 = ps_data.midpoint[1] + 0.5 * ps_data.ds[1]
    y1 = ps_data.midpoint[2] + 0.5 * ps_data.ds[2]
    y2 = ps_data.midpoint[2] + 0.5 * ps_data.ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_mesh_points!(face::Face{Nothing}, x::AV, y::AV)
    ps_data = face.data
    collect_mesh_points_kernal!(Val(Int(face.faceid)), ps_data, x, y)
end
function collect_mesh_points!(::Face, ::AV, ::AV)
    return nothing
end
function collect_save_mesh(DVM_data::DVM_Data)
    x = Vector{Float64}(undef, 0)
    y = Vector{Float64}(undef, 0)
    for i in eachindex(DVM_data.faces)
        face = DVM_data.faces[i]
        collect_mesh_points!(face, x, y)
    end
    return Mesh_Points(x, y)
end
function collect_mesh_points(DVM_data::DVM_Data)
    x = Vector{Float64}(undef, 0)
    y = Vector{Float64}(undef, 0)
    for i in eachindex(DVM_data.faces)
        face = DVM_data.faces[i]
        collect_mesh_points!(face, x, y)
    end
    # return Mesh_Points(
    #             NTuple{length(x),Float64}(x),
    #             NTuple{length(y),Float64}(y)
    #             )
    return Mesh_Points(SVector{length(x),Float64}(x), SVector{length(y),Float64}(y)) # Not a good idea to use SVectors. Change to Point2Point communication later.
end
# function mesh_plot_kernal!(::Val{1},ps_data::AbstractPsData,ax::Axis)
#     x = Vector{Float64}(undef,2);y = Vector{Float64}(undef,2)
#     x[1] = ps_data.midpoint[1] - 0.5*ps_data.ds[1]
#     x[2] = ps_data.midpoint[1] - 0.5*ps_data.ds[1]
#     y[1] = ps_data.midpoint[2] - 0.5*ps_data.ds[2]
#     y[2] = ps_data.midpoint[2] + 0.5*ps_data.ds[2]
#     lines!(ax,x,y,color = :red, linewidth = 0.5)
# end
# function mesh_plot_kernal!(::Val{2},ps_data::AbstractPsData,ax::Axis)
#     x = Vector{Float64}(undef,2);y = Vector{Float64}(undef,2)
#     x[1] = ps_data.midpoint[1] + 0.5*ps_data.ds[1]
#     x[2] = ps_data.midpoint[1] + 0.5*ps_data.ds[1]
#     y[1] = ps_data.midpoint[2] - 0.5*ps_data.ds[2]
#     y[2] = ps_data.midpoint[2] + 0.5*ps_data.ds[2]
#     lines!(ax,x,y,color = :red, linewidth = 0.5)
# end
# function mesh_plot_kernal!(::Val{3},ps_data::AbstractPsData,ax::Axis)
#     x = Vector{Float64}(undef,2);y = Vector{Float64}(undef,2)
#     x[1] = ps_data.midpoint[1] - 0.5*ps_data.ds[1]
#     x[2] = ps_data.midpoint[1] + 0.5*ps_data.ds[1]
#     y[1] = ps_data.midpoint[2] - 0.5*ps_data.ds[2]
#     y[2] = ps_data.midpoint[2] - 0.5*ps_data.ds[2]
#     lines!(ax,x,y,color = :red, linewidth = 0.5)
# end
# function mesh_plot_kernal!(::Val{4},ps_data::AbstractPsData,ax::Axis)
#     x = Vector{Float64}(undef,2);y = Vector{Float64}(undef,2)
#     x[1] = ps_data.midpoint[1] - 0.5*ps_data.ds[1]
#     x[2] = ps_data.midpoint[1] + 0.5*ps_data.ds[1]
#     y[1] = ps_data.midpoint[2] + 0.5*ps_data.ds[2]
#     y[2] = ps_data.midpoint[2] + 0.5*ps_data.ds[2]
#     lines!(ax,x,y,color = :red, linewidth = 0.5)
# end
# function mesh_plot(DVM_data::DVM_Data)
#     f = Figure()
#     ax = Axis(f[1, 1],xlabel = L"x",ylabel = L"y")
#     for i in eachindex(DVM_data.faces)
#         face = DVM_data.faces[i]
#         mesh_plot_kernal!(Val(Int(face.faceid)),face.data,ax)
#     end
#     save("mesh.png",f)
# end

function collect_save_results(p4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    pp = PointerWrapper(p4est)
    field_data = Vector{Field}(undef, pp.local_num_quadrants[])
    index = 1
    for i in eachindex(DVM_data.trees.data)
        for j in eachindex(DVM_data.trees.data[i])
            ps_data = DVM_data.trees.data[i][j]
            field_data[index] = Field(ps_data.midpoint, ps_data.w, ps_data.prim, ps_data.qf)
            index += 1
        end
    end
    return field_data
end
function write_result!(p4est::Ptr{p4est_t}, DVM_data::DVM_Data, name::String)
    results = collect_save_results(p4est, DVM_data)
    JLD2.save_object(name * "_" * string(MPI.Comm_rank(MPI.COMM_WORLD)) * ".jld", results)
end

function read_results(name::String, comm_size::Int)
    results = JLD2.load_object(name * "_0.jld")
    for i = 1:comm_size-1
        append!(results, JLD2.load_object(name * "_" * string(i) * ".jld"))
    end
    return results
end

# function interpolate_result(results::Vector,Var::Symbol)
#     points = Matrix{Float64}(undef,length(results),DIM)
#     values = Vector{Float64}(undef,length(results))
#     for i in eachindex(results)
#         points[i,:] = results[i].midpoint
#         values[i] = getfield(results[i],Var)
#     end
#     scipy.interpolate.CloughTocher2DInterpolator(points,values)
# end

function collect_vs_mesh_points_kernal!(::Val{1}, midpoint::AV, ds::AV, x::AV, y::AV)
    x1 = midpoint[1] - 0.5 * ds[1]
    x2 = midpoint[1] - 0.5 * ds[1]
    y1 = midpoint[2] - 0.5 * ds[2]
    y2 = midpoint[2] + 0.5 * ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_vs_mesh_points_kernal!(::Val{2}, midpoint::AV, ds::AV, x::AV, y::AV)
    x1 = midpoint[1] + 0.5 * ds[1]
    x2 = midpoint[1] + 0.5 * ds[1]
    y1 = midpoint[2] - 0.5 * ds[2]
    y2 = midpoint[2] + 0.5 * ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_vs_mesh_points_kernal!(::Val{3}, midpoint::AV, ds::AV, x::AV, y::AV)
    x1 = midpoint[1] - 0.5 * ds[1]
    x2 = midpoint[1] + 0.5 * ds[1]
    y1 = midpoint[2] - 0.5 * ds[2]
    y2 = midpoint[2] - 0.5 * ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_vs_mesh_points_kernal!(::Val{4}, midpoint::AV, ds::AV, x::AV, y::AV)
    x1 = midpoint[1] - 0.5 * ds[1]
    x2 = midpoint[1] + 0.5 * ds[1]
    y1 = midpoint[2] + 0.5 * ds[2]
    y2 = midpoint[2] + 0.5 * ds[2]
    push!(x, x1, x2)
    push!(y, y1, y2)
end
function collect_vs_mesh_points!(midpoint::AV, ds::AV, x::AV, y::AV)
    for i = 1:4
        collect_vs_mesh_points_kernal!(Val(i), midpoint, ds, x, y)
    end
end
function collect_vs_mesh(ds::AV, midpoint::AM, level::AV)
    x = Vector{Float64}(undef, 0)
    y = Vector{Float64}(undef, 0)
    for i in axes(midpoint, 1)
        collect_vs_mesh_points!(@view(midpoint[i, :]), ds ./ 2^level[i], x, y)
    end
    return x, y
end
function save_result(DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    data_set = DVM_data.data_set
    global_data.gas.sim_time >= data_set.PS_time_interval * global_data.ps_save_index &&
        PS_save(DVM_data)
    global_data.gas.sim_time >= data_set.VS_end_time && return
    global_data.gas.sim_time >= data_set.VS_time_interval * global_data.vs_save_index &&
        VS_save(DVM_data)
end
function save_set(DVM_data::DVM_Data)
    dir_path = "./result/"
    !isdir(dir_path) && mkpath(dir_path)
    data_set = DVM_data.data_set
    result_set = Result_Set(
        DVM_data.global_data,
        data_set.PS_time_interval,
        data_set.VS_time_interval,
        data_set.end_time,
        data_set.VS_end_time,
        MPI.Comm_size(MPI.COMM_WORLD),
    )
    JLD2.save_object(dir_path * "result_set.jld", result_set)
end
function PS_save(DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    field = collect_save_results(DVM_data.p4est, DVM_data)
    mesh = collect_save_mesh(DVM_data)
    dir_path = "./result/PS_result/" * string(global_data.ps_save_index)
    !isdir(dir_path) && mkpath(dir_path)
    ps_result = PS_Result(global_data.gas.sim_time, field, mesh)
    JLD2.save_object(
        dir_path * "/" * string(MPI.Comm_rank(MPI.COMM_WORLD)) * ".jld",
        ps_result,
    )
    global_data.ps_save_index += 1
end
function VS_save(DVM_data::DVM_Data)
    MPI.Comm_rank(MPI.COMM_WORLD) != MPI.Comm_size(MPI.COMM_WORLD) - 1 && return
    global_data = DVM_data.global_data
    ps_data = DVM_data.trees.data[end][end]
    vs_data = ps_data.vs_data
    dir_path = "./result/VS_result/" * string(global_data.vs_save_index)
    !isdir(dir_path) && mkpath(dir_path)
    vs_result = VS_Result(
        global_data.gas.sim_time,
        ps_data.midpoint,
        ps_data.ds,
        vs_data.level,
        vs_data.midpoint,
        vs_data.df,
    )
    JLD2.save_object(
        dir_path * "/" * string(MPI.Comm_rank(MPI.COMM_WORLD)) * ".jld",
        vs_result,
    )
    global_data.vs_save_index += 1
end
function save_VS_final(DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    dir_path = "./result/VS_result/end"
    !isdir(dir_path) && mkpath(dir_path)
    VS_results = Vector{VS_Result}(undef, 0)
    for i in eachindex(DVM_data.trees.data)
        for j in eachindex(DVM_data.trees.data[i])
            ps_data = DVM_data.trees.data[i][j]
            if ps_data.midpoint[2] < 0.7 + ps_data.ds[2]
                vs_data = ps_data.vs_data
                vs_result = VS_Result(
                    global_data.gas.sim_time,
                    ps_data.midpoint,
                    ps_data.ds,
                    vs_data.level,
                    vs_data.midpoint,
                    vs_data.df,
                )
                push!(VS_results, vs_result)
            end
        end
    end
    JLD2.save_object(
        dir_path * "/" * string(MPI.Comm_rank(MPI.COMM_WORLD)) * ".jld",
        VS_results,
    )
end
function result2animation()
    result_set = JLD2.load_object("./result/result_set.jld")
    global_data = result_set.global_data
    mpi_size = result_set.mpi_size
    nframes = Int(round(result_set.end_time / result_set.PS_time_interval))
    framerate = Int(round(1 / result_set.PS_time_interval))
    x, y, mesh_xs, mesh_ys, variables = read_result(global_data, nframes, mpi_size)
    index = Observable(0)
    f = Figure()
    ax = Axis(f[1, 1])
    variable = @lift(variables[$index+1])
    co = contourf!(x, y, variable)
    Colorbar(f[1, 2], co)
    @lift begin
        mesh_x = mesh_xs[$index+1]
        mesh_y = mesh_ys[$index+1]
        mesh_plot(mesh_x, mesh_y, ax)
    end
    index_iter = 0:nframes

    record(f, "PS_field.mp4", index_iter; framerate = framerate) do i
        # x,y,variable = read_resulte_at_index(i,DVM_data.global_data)
        # Axis(f[1, 1])
        # co = contourf!(x,y,variable)
        # Colorbar(f[1, 2], co)
        index[] = i
        @show i
    end
end
function read_result_at_index!(
    index::Int,
    global_data::Global_Data,
    mesh_xs::Vector{Vector{Float64}},
    mesh_ys::Vector{Vector{Float64}},
    variables::Vector{Matrix{Float64}},
    mpi_size::Int,
)
    fields = Vector{Field}(undef, 0)
    mesh_x = Vector{Float64}(undef, 0)
    mesh_y = Vector{Float64}(undef, 0)
    for j = 0:mpi_size-1
        PS_results = JLD2.load_object(
            "./result/PS_result/" * string(index) * "/" * string(j) * ".jld",
        )
        append!(fields, PS_results.field)
        append!(mesh_x, PS_results.mesh.x)
        append!(mesh_y, PS_results.mesh.y)
    end
    x, y, variable = reshape_solutions(fields, global_data, :prim, 4)
    variable .= 1 ./ variable
    push!(variables, variable)
    push!(mesh_xs, mesh_x)
    push!(mesh_ys, mesh_y)
    return x, y, variable
end
function read_result(global_data::Global_Data, nframes::Int, mpi_size::Int)
    mesh_xs = Vector{Vector{Float64}}(undef, 0)
    mesh_ys = Vector{Vector{Float64}}(undef, 0)
    variables = Vector{Matrix{Float64}}(undef, 0)
    x, y, _ = read_result_at_index!(0, global_data, mesh_xs, mesh_ys, variables, mpi_size)
    for i = 1:nframes
        read_result_at_index!(i, global_data, mesh_xs, mesh_ys, variables, mpi_size)
    end
    return x, y, mesh_xs, mesh_ys, variables
end
# function pre_read_resulte_at_index(index::Int,global_data::Global_Data)
#     fields = Vector{Field}(undef,0)
#     for j = 0:5
#         PS_results = JLD2.load_object("./result/PS_result/"*string(index)*"/"*string(j)*".jld")
#         append!(fields,PS_results.field)
#     end
#     x, y, _ = reshape_solutions(fields,global_data,:prim,4)
#     return x,y
# end
function reshape_solutions_vs(midpoint::AM, df::AM, global_data::Global_Data)
    Nx, Ny = global_data.vs_trees_num
    xmin, xmax, ymin, ymax = global_data.quadrature
    z = @view(df[:, 1])
    itp = scipy.interpolate.CloughTocher2DInterpolator(midpoint, PyArray(z))
    dx = (xmax - xmin) / Nx / 2^DVM_VS_MAXLEVEL
    dy = (ymax - ymin) / Ny / 2^DVM_VS_MAXLEVEL
    X = collect(xmin+EPS+dx/2:dx:xmax-EPS-dx/2)
    Y = collect(ymin+EPS+dy/2:dy:ymax-EPS-dy/2)
    pyZ = np.zeros((length(X), length(Y)))
    for i in eachindex(X)
        for j in eachindex(Y)
            pyZ[i-1, j-1] = itp(X[i], Y[j])[]
        end
    end
    Z = pyconvert(Matrix{Float64}, pyZ)
    return X, Y, Z
end
function read_vs_result(global_data::Global_Data, nframes::Int, mpi_size::Int)
    mesh_xs = Vector{Vector{Float64}}(undef, 0)
    mesh_ys = Vector{Vector{Float64}}(undef, 0)
    hs = Vector{Matrix{Float64}}(undef, 0)
    x, y, _ = read_vs_result_at_index!(0, global_data, mesh_xs, mesh_ys, hs, mpi_size)
    for i = 1:nframes-1
        read_vs_result_at_index!(i, global_data, mesh_xs, mesh_ys, hs, mpi_size)
    end
    return x, y, mesh_xs, mesh_ys, hs
end
function read_vs_result_at_index!(
    index::Int,
    global_data::Global_Data,
    mesh_xs::Vector{Vector{Float64}},
    mesh_ys::Vector{Vector{Float64}},
    hs::Vector{Matrix{Float64}},
    mpi_size::Int,
)
    VS_results = JLD2.load_object(
        "./result/VS_result/" * string(index) * "/" * string(mpi_size - 1) * ".jld",
    )
    ds = zeros(DIM)
    for j = 1:DIM
        ds[j] =
            (global_data.quadrature[2*j] - global_data.quadrature[2*j-1]) /
            global_data.vs_trees_num[j]
    end
    x, y, h = reshape_solutions_vs(VS_results.midpoint, VS_results.df, global_data)
    mesh_x, mesh_y = collect_vs_mesh(ds, VS_results.midpoint, VS_results.level)
    push!(hs, h)
    push!(mesh_xs, mesh_x)
    push!(mesh_ys, mesh_y)
    return x, y, h
end
function vsresult2animation()
    result_set = JLD2.load_object("./result/result_set.jld")
    global_data = result_set.global_data
    mpi_size = result_set.mpi_size
    nframes = Int(round(result_set.VS_end_time / result_set.VS_time_interval))
    framerate = Int(round(1 / result_set.VS_time_interval) / 2)
    x, y, mesh_xs, mesh_ys, hs = read_vs_result(global_data, nframes, mpi_size)
    index = Observable(0)
    f = Figure()
    ax = Axis(f[1, 1])
    h = @lift(hs[$index+1])
    co = contourf!(x, y, h)
    Colorbar(f[1, 2], co)
    @lift begin
        mesh_x = mesh_xs[$index+1]
        mesh_y = mesh_ys[$index+1]
        mesh_plot(mesh_x, mesh_y, ax)
    end
    index_iter = 0:nframes-1

    record(f, "VS_field.mp4", index_iter; framerate = framerate) do i
        index[] = i
        @show i
    end
end
function mat_3d()
    result_set = JLD2.load_object("./result/result_set.jld")
    global_data = result_set.global_data
    ds = zeros(DIM)
    for j = 1:DIM
        ds[j] =
            (global_data.quadrature[2*j] - global_data.quadrature[2*j-1]) /
            global_data.vs_trees_num[j]
    end
    mpi_size = result_set.mpi_size
    u = Vector{Float64}(undef, 0)
    v = Vector{Float64}(undef, 0)
    x = Vector{Float64}(undef, 0)
    h = Vector{Float64}(undef, 0)
    mesh_u = Vector{Float64}(undef, 0)
    mesh_v = Vector{Float64}(undef, 0)
    mesh_x = Vector{Float64}(undef, 0)
    dir_path = "./result/VS_result/end"
    for i = 0:mpi_size-1
        vs_results = JLD2.load_object(dir_path * "/" * string(i) * ".jld")
        isempty(vs_results) && continue
        for j in eachindex(vs_results)
            append!(h, @view(vs_results[j].df[:, 1]))
            append!(u, @view(vs_results[j].midpoint[:, 1]))
            append!(v, @view(vs_results[j].midpoint[:, 2]))
            append!(x, vs_results[j].ps_midpoint[1] * ones(size(vs_results[j].midpoint, 1)))
            collect_mesh_3d!(vs_results[j], ds, mesh_u, mesh_v, mesh_x)
        end
    end
    Nx, _ = global_data.trees_num
    xmin, xmax, _, _ = global_data.geometry
    dx = (xmax - xmin) / Nx / 2^DVM_PS_MAXLEVEL
    gridX = collect(xmin+dx/2:dx:xmax-dx/2)
    Nu, Nv = global_data.vs_trees_num
    umin, umax, vmin, vmax = global_data.quadrature
    du = (umax - umin) / Nu / 2^DVM_VS_MAXLEVEL
    dv = (vmax - vmin) / Nv / 2^DVM_VS_MAXLEVEL
    gridU = collect(umin+du/2:du:umax-du/2)
    gridV = collect(vmin+dv/2:dv:vmax-dv/2)
    file = matopen("mat_3d.mat", "w")
    write(file, "u", u)
    write(file, "v", v)
    write(file, "x", x)
    write(file, "h", h)
    write(file, "gridX", gridX)
    write(file, "gridU", gridU)
    write(file, "gridV", gridV)
    write(file, "mesh_x", mesh_x)
    write(file, "mesh_u", mesh_u)
    write(file, "mesh_v", mesh_v)
    close(file)
end
function collect_mesh_3d!(vs_result, ds, mesh_x, mesh_y, mesh_z)
    midpoint = vs_result.midpoint
    for j in axes(midpoint, 1)
        for i = 1:12
            collect_mesh_3d!(
                Val(i),
                midpoint[j, 1],
                midpoint[j, 2],
                vs_result.ps_midpoint[1],
                ds[1] / 2^vs_result.level[j],
                ds[2] / 2^vs_result.level[j],
                vs_result.ps_ds[1],
                mesh_x,
                mesh_y,
                mesh_z,
            )
        end
    end
end
function collect_mesh_3d!(::Val{1}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y - 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z - 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{2}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x + 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z - 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{3}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y + 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z - 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{4}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x - 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z - 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{5}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y - 0.5 * dy
    z1 = z + 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{6}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x + 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z + 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{7}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y + 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z + 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{8}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x - 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z + 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{9}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x - 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y - 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{10}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x + 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y - 0.5 * dy
    y2 = y - 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{11}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x + 0.5 * dx
    x2 = x + 0.5 * dx
    y1 = y + 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end
function collect_mesh_3d!(::Val{12}, x::Float64, y, z, dx, dy, dz, X::AV, Y::AV, Z::AV)
    x1 = x - 0.5 * dx
    x2 = x - 0.5 * dx
    y1 = y + 0.5 * dy
    y2 = y + 0.5 * dy
    z1 = z - 0.5 * dz
    z2 = z + 0.5 * dz
    push!(X, x1, x2)
    push!(Y, y1, y2)
    push!(Z, z1, z2)
end

function plot_at_index(i::Int)
    result_set = JLD2.load_object("./result/result_set.jld")
    global_data = result_set.global_data
    mpi_size = result_set.mpi_size
    mesh_xs = Vector{Vector{Float64}}(undef, 0)
    mesh_ys = Vector{Vector{Float64}}(undef, 0)
    variables = Vector{Matrix{Float64}}(undef, 0)
    x, y, _ = read_result_at_index!(i, global_data, mesh_xs, mesh_ys, variables, mpi_size)
    f = Figure()
    ax = Axis(f[1, 1])
    co = contourf!(x, y, variables[1])
    Colorbar(f[1, 2], co)
    save("PS_field.png", f)

    f = Figure()
    ax = Axis(f[1, 1])
    mesh_plot(mesh_xs[1], mesh_ys[1], ax)
    save("mesh.png", f)
end
