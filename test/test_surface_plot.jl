using KitAMR
boundary_solutions = KitAMR.Boundary_Solution(Vector{Float64}[],Vector{Float64}[],KitAMR.Boundary_PS_Solution[])
path = "./result2026-03-02_03-23/"
for i in 1:180
    bs = KitAMR.load_object(path*"boundary_result_"*string(i-1)*".jld2")
    append!(boundary_solutions.midpoints,bs[1].midpoints)
    append!(boundary_solutions.normal,bs[1].normal)
    append!(boundary_solutions.ps_solutions,bs[1].ps_solutions)
end
points = [boundary_solutions.midpoints[i][j] for i in eachindex(boundary_solutions.midpoints), j in 1:3] |> permutedims
cells = [KitAMR.MeshCell(KitAMR.VTKCellTypes.VTK_VERTEX,[i]) for i in eachindex(boundary_solutions.midpoints)]
KitAMR.vtk_grid("surf_test",points,cells) do vtk
    vtk["rho"] = [x.prim[1] for x in boundary_solutions.ps_solutions]
    vtk["U"] = [x.prim[2] for x in boundary_solutions.ps_solutions]
    vtk["V"] = [x.prim[3] for x in boundary_solutions.ps_solutions]
    vtk["W"] = [x.prim[4] for x in boundary_solutions.ps_solutions]
    vtk["T"] = [1.0/x.prim[5] for x in boundary_solutions.ps_solutions]
    vtk["p11"] = [x.p[1] for x in boundary_solutions.ps_solutions]
    vtk["p12"] = [x.p[2] for x in boundary_solutions.ps_solutions]
    vtk["p13"] = [x.p[3] for x in boundary_solutions.ps_solutions]
    vtk["p22"] = [x.p[4] for x in boundary_solutions.ps_solutions]
    vtk["p23"] = [x.p[5] for x in boundary_solutions.ps_solutions]
    vtk["p33"] = [x.p[6] for x in boundary_solutions.ps_solutions ]
    vtk["normal"] = ([x[1] for x in boundary_solutions.normal],[x[2] for x in boundary_solutions.normal],
        [x[3] for x in boundary_solutions.normal])
end