using CairoMakie, JLD2

cd(@__DIR__)
include("DVM.jl")

result_set = JLD2.load_object("result/PS_result/8/1.jld")

using JLD
JLD.load("result/PS_result/8/1.jld")



global_data = result_set.global_data
mpi_size = result_set.mpi_size
mesh_xs = Vector{Vector{Float64}}(undef, 0)
mesh_ys = Vector{Vector{Float64}}(undef, 0)
variables = Vector{Matrix{Float64}}(undef, 0)
i = 1
x, y, _ = read_result_at_index!(i, global_data, mesh_xs, mesh_ys, variables, mpi_size)
f = Figure()
ax = Axis(f[1, 1])
co = contourf!(x, y, variables[1])
Colorbar(f[1, 2], co)
save("test_vs_adaptive.png", f)


result_set = JLD2.load_object("./result/result_set.jld")
global_data = result_set.global_data
mpi_size = result_set.mpi_size
mesh_xs = Vector{Vector{Float64}}(undef, 0)
mesh_ys = Vector{Vector{Float64}}(undef, 0)
variables = Vector{Matrix{Float64}}(undef, 0)
i = 1
x, y, _ = read_result_at_index!(i, global_data, mesh_xs, mesh_ys, variables, mpi_size)
f = Figure()
ax = Axis(f[1, 1])
co = contourf!(x, y, variables[1])
Colorbar(f[1, 2], co)
save("test_vs_adaptive.png", f)

f = Figure()
ax = Axis(f[1, 1])
mesh_plot(mesh_xs[1], mesh_ys[1], ax)
save("mesh.png", f)



using Plots, JLD



d = JLD.load("result/PS_result/8/1.jld")

ps = d["single_stored_object"]
