using KitAMR,MPI

MPI.Init()
config = KitAMR.read_config("configure_2D.txt")
ps4est,amr = KitAMR.init(config);

# fp = PointerWrapper(ps4est)
# ps_solution = Vector{KitAMR.PS_Solution}(undef,fp.local_num_quadrants[])
# neighbor_nums = Vector{Vector{Int}}(undef,fp.local_num_quadrants[])
# trees = amr.field.trees.data
# config = amr.global_data.config
# index = 1
# for i in eachindex(trees)
#     for j in eachindex(trees[i])
#         ps_data = trees[i][j]
#         ps_solution[index] = KitAMR.PS_Solution(ps_data)
#         neighbor_nums[index] = KitAMR.neighbor_num(ps_data,ps4est,amr,index-1)
#         index+=1
#     end
# end
# solution = KitAMR.Solution(ps_solution)
# rank = MPI.Comm_rank(MPI.COMM_WORLD)
# result = KitAMR.Result(solution,KitAMR.MeshInfo(neighbor_nums))
# dir_path = "./result"*string(now())*"/"
# !isdir(dir_path) && mkpath(dir_path)
# p4est_save("p",ps4est,Cint(0))
# if rank==0
#     size = MPI.Comm_size(MPI.COMM_WORLD)
#     solverset = KitAMR.SolverSet(config,size)
#     save_object(dir_path * "solverset.jld2", solverset)
# end
# save_object(dir_path * "result_"*string(rank)*".jld2", result)



# cnn = Ptr{Ptr{KitAMR.p4est_connectivity_t}}(Libc.malloc(sizeof(Ptr{Ptr{p4est_connectivity_t}})))
# ps4est_reload = p4est_load("p",MPI.COMM_WORLD,Cint(0),Cint(0),C_NULL,cnn)
# path = "./"*dirname
# solverset = load_object(path*"/solverset.jld2")
# result = nothing
# for i = 1:solverset.mpi_size
#     if i==1
#         result = load_object(path*"/result_"*string(i-1)*".jld2")
#     else
#         r = load_object(path*"/result_"*string(i-1)*".jld2")
#         for j in (solution.ps_data,solution.vs_data,mesh_info.neighbor_nums)
#             @eval append!($result.$j,$r.$j)
#         end
#     end
# end

# vtk_cnn = Vector{Vector{Int}}(undef,length(result.solution.ps_solutions))
# neighbor_nums = result.mesh_info.neighbor_nums
# for i in eachindex(neighbor_nums)
#     for j in eachindex(neighbor_nums[i])
#         neighbor_nums[i][j]==0&&(neighbor_nums[i][j]=1)
#     end
# end
# for i in eachindex(vtk_cnn)
#     addi_points_num = 0
#     for j in eachindex(neighbor_nums[i])
#         addi_points_num+= neighbor_nums[i][j]==2^(DIM-1) ? 1 : 0
#     end
#     vtk_cnn[i] = Vector{Int}(undef,2^DIM+addi_points_num)
# end
# points = Vector{Float64}[]
# data = Vector{Any}()
# for el in (points,vtk_cnn,neighbor_nums)
#     push!(data,el)
# end
# p_data = pointer_from_objref(data)
# GC.@preserve data AMR_corner_iterate(ps4est;user_data = p_data) do ip,data
#     points,vtk_cnn,neighbor_nums = unsafe_pointer_to_objref(data)
#     DIM=isa(ip,PointerWrapper{p4est_iter_corner_info_t}) ? 2 : 3
#     for i in 1:ip.sides.elem_count[]
#         side = KitAMR.iPointerWrapper(ip.sides,KitAMR.p4est_iter_corner_side_t,i-1)
#         cornerid = side.corner[]+1 # z-order
#         if i==1
#             ds,midpoint = KitAMR.quad_to_cell(ip.p4est,side.treeid[],side.quad)
#             point = @. midpoint+0.5*ds*KitAMR.RMT[DIM][cornerid]
#             push!(points,point)
#         end
#         id = length(points)
#         quadid = KitAMR.local_quadid(ip,side)
#         neighbor_num = neighbor_nums[quadid+1]
#         cornerid==1&&(vtk_cnn[quadid+1][1]=id)
#         cornerid==2&&(vtk_cnn[quadid+1][neighbor_num[3]+1]=id)
#         cornerid==3&&(vtk_cnn[quadid+1][sum(neighbor_num[2:4])+1]=id)
#         cornerid==4&&(vtk_cnn[quadid+1][sum(neighbor_num[2:3])+1]=id)
#     end
#     return nothing
# end
# GC.@preserve data AMR_face_iterate(ps4est;user_data = p_data) do ip,data
#     points,vtk_cnn,neighbor_nums = unsafe_pointer_to_objref(data)
#     DIM=isa(ip,PointerWrapper{p4est_iter_face_info_t}) ? 2 : 3
#     ip.sides.elem_count[]==1&&return nothing
#     side1 = KitAMR.iPointerWrapper(ip.sides,KitAMR.p4est_iter_face_side_t,0)
#     side2 = KitAMR.iPointerWrapper(ip.sides,KitAMR.p4est_iter_face_side_t,1)
#     (side1.is_hanging[]==0&&side2.is_hanging[]==0)&&return nothing
#     for side in (side1,side2)
#             if side.is_hanging[]==0
#                 ds,midpoint = KitAMR.quad_to_cell(ip.p4est,side.treeid[],side.is.full.quad)
#                 faceid = side.face[]+1
#                 point = @. midpoint+0.5*ds*KitAMR.NMT[DIM][faceid]
#                 push!(points,point)
#                 id = length(points)
#                 quadid = KitAMR.local_quadid(ip.p4est,side.treeid[],side.is.full.quadid[])
#                 neighbor_num = neighbor_nums[quadid+1]
#                 faceid==1&&(@inbounds vtk_cnn[quadid+1][end]=id)
#                 faceid==2&&(@inbounds vtk_cnn[quadid+1][neighbor_num[3]+2]=id)
#                 faceid==3&&(@inbounds vtk_cnn[quadid+1][2]=id)
#                 faceid==4&&(@inbounds vtk_cnn[quadid+1][sum(neighbor_num[2:3])+2]=id)
#             end
#     end
#     for side in (side1,side2)
#             if side.is_hanging[]==1
#                 faceid = side.face[]+1
#                 id = length(points)
#                 quadids = side.is.hanging.quadid[]
#                 quadid1 = KitAMR.local_quadid(ip.p4est,side.treeid[],quadids[1])
#                 quadid2 = KitAMR.local_quadid(ip.p4est,side.treeid[],quadids[2])
#                 neighbor_num1 = neighbor_nums[quadid1+1]
#                 neighbor_num2 = neighbor_nums[quadid2+1]
#                 faceid==1&&(@inbounds vtk_cnn[quadid1+1][end]=id;vtk_cnn[quadid2+1][1]=id)
#                 faceid==2&&(@inbounds vtk_cnn[quadid1+1][sum(neighbor_num1[2:3])+1]=id;vtk_cnn[quadid2+1][neighbor_num2[3]+1]=id)
#                 faceid==3&&(@inbounds vtk_cnn[quadid1+1][2]=id;vtk_cnn[quadid2+1][1]=id)
#                 faceid==4&&(@inbounds vtk_cnn[quadid1+1][sum(neighbor_num1[2:3])+1]=id;vtk_cnn[quadid2+1][end]=id)
#             end
#     end
#     return nothing
# end
# cells = [MeshCell(PolyData.Polys(),cnn) for cnn in vtk_cnn]
# vertices = Matrix{Float64}(undef,2,length(points))
# for i in eachindex(points)
#     @inbounds vertices[:,i] .= points[i]
# end
# vtk_grid("test_vtk_output",vertices,cells) do vtk
#     vtk["T"] = [1/ps_solution.prim[end] for ps_solution in result.solution.ps_solutions]
# end

# GC.@preserve data AMR_face_iterate(ps4est;user_data = p_data) do ip,data
#     points,vtk_cnn,neighbor_nums = unsafe_pointer_to_objref(data)
#     DIM=isa(ip,PointerWrapper{p4est_iter_face_info_t}) ? 2 : 3
#     ip.sides.elem_count[]==1&&return nothing
#     side1 = KitAMR.iPointerWrapper(ip.sides,KitAMR.p4est_iter_face_side_t,0)
#     side2 = KitAMR.iPointerWrapper(ip.sides,KitAMR.p4est_iter_face_side_t,1)
#     (side1.is_hanging[]==0&&side2.is_hanging[]==0)&&return nothing
#     for side in (side1,side2)
#         @eval begin
#             if $side.is_hanging[]==0
#                 ds,midpoint = KitAMR.quad_to_cell($ip.p4est,$side.treeid[],$side.is.full.quad)
#                 faceid = $side.face[]+1
#                 point = @. midpoint+0.5*ds*KitAMR.NMT[DIM][faceid]
#                 push!($points,point)
#                 id = length($points)
#                 quadid = KitAMR.local_quadid($ip.p4est,$side.treeid[],$side.is.full.quadid[])
#                 neighbor_num = $neighbor_nums[quadid+1]
#                 faceid==1&&(@inbounds $vtk_cnn[quadid+1][end]=id)
#                 faceid==2&&(@inbounds $vtk_cnn[quadid+1][neighbor_num[3]+2]=id)
#                 faceid==3&&(@inbounds $vtk_cnn[quadid+1][2]=id)
#                 faceid==4&&(@inbounds $vtk_cnn[quadid+1][sum(neighbor_num[2:3])+2]=id)
#             end
#         end
#     end
#     for side in (side1,side2)
#         @eval begin
#             if $side.is_hanging[]==1
#                 faceid = $side.face[]+1
#                 id = length(points)
#                 quadids = $side.is.hanging.quadid[]
#                 quadid1 = KitAMR.local_quadid($ip.p4est,$side.treeid[],quadids[1])
#                 quadid2 = KitAMR.local_quadid($ip.p4est,$side.treeid[],quadids[2])
#                 neighbor_num1 = $neighbor_nums[quadid1+1]
#                 neighbor_num2 = $neighbor_nums[quadid2+1]
#                 faceid==1&&(@inbounds $vtk_cnn[quadid1+1][end]=id;$vtk_cnn[quadid2+1][1]=id)
#                 faceid==2&&(@inbounds $vtk_cnn[quadid1+1][sum(neighbor_num1[2:3])+1]=id;$vtk_cnn[quadid2+1][neighbor_num2[3]+1]=id)
#                 faceid==3&&(@inbounds $vtk_cnn[quadid1+1][2]=id;$vtk_cnn[quadid2+1][1]=id)
#                 faceid==4&&(@inbounds $vtk_cnn[quadid1+1][sum(neighbor_num1[2:3])+1]=id;$vtk_cnn[quadid2+1][end]=id)
#             end
#         end
#     end
#     return nothing
# end
