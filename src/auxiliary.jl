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
function gaussian_area(A::AbstractMatrix) # DIMxN
    area = 0
    for i in axes(A,2)
        j = i%size(A,2)+1
        area+=@views det(A[:,[i,j]])
    end
    return 0.5*abs(area)
end
function cut_rect(n::Vector{Float64},vertices::Vector{Vector{Float64}}) # The splitting line at the boundary with normal direction n. Return the areas of the two part. The former one locates at the direction n.
    # The cases where n is oriented to along the axes are excluded in advance.
    #=
    y
#     7(8)|-----6-----|5(6)
#         |     |     |
#         8-----|-----4
#         |     |     |
#     1(2)|-----2-----|3(4) x
    =#
    points = Vector{Float64}[];indices = Int[];i = 0
    xmin,xmax,ymin,ymax = vertices[1][1],vertices[2][1],vertices[1][2],vertices[3][2]
    y = -n[1]/n[2]*xmin
    if y<ymax&&y>ymin
        push!(points,[xmin,y])
        push!(indices,8)
        i+=1
    end
    y = -n[1]/n[2]*xmax
    if y<ymax&&y>ymin
        push!(points,[xmax,y])
        push!(indices,4)
        i+=1
    end
    x = -n[2]/n[1]*ymin
    if x<xmax&&x>xmin
        push!(points,[x,ymin])
        push!(indices,2)
        i+=1
    end
    x = -n[2]/n[1]*ymax
    if x<xmax&&x>xmin
        push!(points,[x,ymax])
        push!(indices,6)
        i+=1
    end
    clp = findall(x->dot(x,n)==0.,vertices) # Number of vertices lying on the discontinuity
    if length(clp)+i<2
        return false,0.,0.
    else
        append!(points,vertices)
        append!(indices,CLP) # CLP is defined in dim.jl
        sid = sortperm(indices)
        indices = indices[sid]
        points = points[sid]
        aid = findfirst(x->iseven(x)||in(x,CLP[clp]),indices)
        bid = findlast(x->iseven(x)||in(x,CLP[clp]),indices)
        A = Matrix{Float64}(undef,2,bid-aid+1);B = Matrix{Float64}(undef,2,length(indices)-(bid-aid)+1)
        for i in aid:bid
            A[:,i-aid+1] .= points[i]
        end
        for i in 0:length(indices)-(bid-aid)
            index = bid+i>length(indices) ? (bid+i)%length(indices) : bid+i
            B[:,i+1] .= points[index]
        end
        l = points[bid]-points[aid]
        if n[1]*l[2]-n[2]*l[1]<0
            gas_weight = gaussian_area(A)
            return true,gas_weight,((xmax-xmin)*(ymax-ymin))-gas_weight # solid first
        else
            solid_weight = gaussian_area(A)
            return true,((xmax-xmin)*(ymax-ymin))-solid_weight,solid_weight
        end
    end
end
function ray_casting(point::Vector{Float64},vertices::Vector{Vector{Float64}}) # Inside the closed curve?
    n = length(vertices)
    count = 0
    for i in 1:n
        j = (i % n) + 1
        x_i, y_i = vertices[i]
        x_j, y_j = vertices[j]
        px, py = point
        if (py > min(y_i, y_j)) && (py <= max(y_i, y_j))
            if y_i == y_j 
                continue
            end
            x_intersect = (py-y_i) * (x_j - x_i) / (y_j - y_i) + x_i
            if x_intersect > px && (x_intersect <= max(x_i, x_j))
                count +=1
            end
        end
    end
    return (count % 2) == 1
end
function calc_normal(p1,p2,s1,s2) # Calculate the unit normal vector of the vector (s2-s1), in the direction of (p2-p1).
    n = [s2[2]-s1[2],s1[1]-s2[1]]
    n /= norm(n)
    if dot(n,p2-p1)<0
        n = -n
    end
    return n
end
function find_horizontal_intersection_2D(s_midpoint, f_midpoint, points)
    n = length(points)
    x_seg_min, x_seg_max = minmax(s_midpoint[1], f_midpoint[1])
    y0 = s_midpoint[2]
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n + 1]
        x1, y1 = p1
        x2, y2 = p2
        if (y1 ≤ y0 ≤ y2) || (y2 ≤ y0 ≤ y1)
            t = (y0 - y1) / (y2 - y1)
            if 0 < t < 1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    return [x_intersect, y0],calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;n0/=norm(n0)
                    return [x_intersect, y0],n0
                end
            end
        end
    end
    @show s_midpoint f_midpoint
    throw("Intersect error!")
    return Float64[],Float64[]
end

function find_horizontal_intersection_3D(s_midpoint, f_midpoint, points)
    n = length(points)
    x_seg_min, x_seg_max = minmax(s_midpoint[1], f_midpoint[1])
    y0 = s_midpoint[2]; z0 = s_midpoint[3]
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n + 1]
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        if ((y1 ≤ y0 ≤ y2) || (y2 ≤ y0 ≤ y1))&&((z1 ≤ z0 ≤ z2) || (z2 ≤ z0 ≤ z1))
            t = (y0 - y1) / (y2 - y1)
            if 0 < t < 1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    return [x_intersect, y0],calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;n0/=norm(n0)
                    return [x_intersect, y0],n0
                end
            end
        end
    end
    @show s_midpoint f_midpoint
    throw("Intersect error!")
    return Float64[],Float64[]
end

function find_vertical_intersection_2D(s_midpoint, f_midpoint, points)
    n = length(points)
    y_seg_min, y_seg_max = minmax(s_midpoint[2], f_midpoint[2])
    x0 = s_midpoint[1]
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n + 1]
        x1, y1 = p1
        x2, y2 = p2
        if (x1 ≤ x0 ≤ x2) || (x2 ≤ x0 ≤ x1)
            t = (x0 - x1) / (x2 - x1)
            if 0 < t < 1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    return [x0, y_intersect],calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;n0/=norm(n0)
                    return [x0, y_intersect],n0
                end
            end
        end
    end
    @show s_midpoint f_midpoint
    throw("Intersect error!")
    return Float64[],Float64[]
end

function find_transverse_intersection_3D(s_midpoint, f_midpoint, points)
    n = length(points)
    z_seg_min, z_seg_max = minmax(s_midpoint[3], f_midpoint[3])
    x0 = s_midpoint[1]; y0 = 
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n + 1]
        x1, y1 = p1
        x2, y2 = p2
        if (x1 ≤ x0 ≤ x2) || (x2 ≤ x0 ≤ x1)
            t = (x0 - x1) / (x2 - x1)
            if 0 < t < 1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    return [x0, y_intersect],calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;n0/=norm(n0)
                    return [x0, y_intersect],n0
                end
            end
        end
    end
    @show s_midpoint f_midpoint
    throw("Intersect error!")
    return Float64[],Float64[]
end

function find_intersections_2D(s_midpoint, f_midpoint, closed_curve)
    x1, y1 = s_midpoint
    x2, y2 = f_midpoint
    if abs(y1-y2)<EPS
        return find_horizontal_intersection(s_midpoint, f_midpoint, closed_curve)
    elseif abs(x1-x2)<EPS
        return find_vertical_intersection(s_midpoint, f_midpoint, closed_curve)
    else
        @show s_midpoint f_midpoint
        throw("Horizontal or vertical line segment is expected!")
    end
end
function find_intersections_3D(s_midpoint, f_midpoint, closed_curve)
    x1, y1, z1 = s_midpoint
    x2, y2, z2 = f_midpoint
    if abs(y1-y2)<EPS && abs(z1-z2)<EPS
        return find_horizontal_intersection(s_midpoint, f_midpoint, closed_curve)
    elseif abs(x1-x2)<EPS && abs(z1-z2)<EPS
        return find_vertical_intersection(s_midpoint, f_midpoint, closed_curve)
    else
        return find_transverse_intersection(s_midpoint,f_midpoint,closed_curve)
    end
end
function fieldvalues_fn(vs_data,aux_df)
    NDF = typeof(vs_data).parameters[2]
    return [aux_df[:,i] for i in 1:NDF]
end
function fieldvalues_fn(vs_data)
    NDF = typeof(vs_data).parameters[2]
    return [vs_data.df[:,i] for i in 1:NDF]
end
function cut_cube_rotate(n::Vector{Float64})
    C = zeros(2,3)
    e1 = [1.,0.,0.]
    e1 = cross(n,e1);e1./=norm(e1)
    e2 = cross(n,e1);e2./=norm(e2)
    if dot(cross(e1,e2),n)>0
        C[1,:].=e1;C[2,:].=e2
    else
        C[2,:].=e1;C[1,:].=e2
    end
    return C
end

"""
Calculate the flux of the vector (x,y,z)/3 through the vertical faces composed by intersecting points and vertices.
Normally, the cases where the connecting line of the two neighboring intersecting pionts are parallel to grid lines should have been got rid of before the function is called.
"""
function Vertical_Volume_Flux(points::AbstractVector{Vector{Float64}},midpoint::Vector{Float64},vertices::AbstractMatrix{Float64})
    H = 0.;N = length(points)
    for i in eachindex(points)
        p1 = points[i%N+1];p2 = points[i]
        dx= p1-p2
        id = findall(x->abs(x)<EPS,dx)
        if length(id)>1
            throw(`Cut-cube Error!`)
        else
            dir = id[1]
            vid = findall(x->abs(x[dir]-points[i][dir])<EPS,eachcol(vertices))
            A = hcat(vertices[:,vid],p1,p2)
            center = vec(mean(A,dims=2))
            id1 = dir%3+1;id2 = (dir+1)%3+1
            phi=@views [atan(x[id2]-center[id2],x[id1]-center[id1]) for x in eachcol(A)]
            id_gauss = sortperm(phi)
            H+=sign(p1[dir]-midpoint[dir])*p1[dir]/3.0*gaussian_area(@views A[[id1,id2],id_gauss])
        end
    end
    return H
end
function cut_cube(n::Vector{Float64},C::Matrix{Float64},midpoint::Vector{Float64},ddu::Vector{Float64},vertices::Matrix{Float64}) # 3.627 μs (215 allocations: 8.45 KiB). Acceptable?
    vltable = [[1,2],[3,4],[7,8],[5,6],[1,3],[2,4],[6,8],[5,7],[1,5],[2,6],[4,8],[3,7]] # vertices-edges table
    points = Vector{Vector{Float64}}(undef,6);index = 1
    dirs = permutedims(vertices)*n
    for i in eachindex(vltable)
        flag = dirs[vltable[i][1]]*dirs[vltable[i][2]] # flag==0: cut any end of the edge; flag<0: cut the edge; flag>0: not cut the edge
        if abs(flag)<EPS
            if cld(i,4)==1 # avoid redundancy
                if abs(dirs[vltable[i][1]])<EPS # end A intersects
                    points[index] = vertices[:,vltable[i][1]];index+=1
                else # end B intersects
                    points[index] = vertices[:,vltable[i][2]];index+=1
                end
            end
        elseif flag<0 # intersects between the two ends
            point = vertices[:,vltable[i][1]]
            dir = cld(i,4);point[dir]=0.
            point[dir] = -dot(point,n)/(n[dir])
            points[index] = point;index+=1
        end
    end
    index<4&&return false,0.,0.
    points = points[1:index-1]
    posid = findall(x->x>EPS,dirs)
    pos = length(posid)<4 ? true : false # H represents solid?
    if any(x->abs(x)<EPS,n) # Simple case
        dir = findfirst(x->abs(x)<EPS,n) 
        pid = findall(x->abs(x[dir]-vertices[dir])<EPS,points) # vertices[dir]: the dir-th component of the first vertex
        if pos
            vid = findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x,n)>EPS,eachcol(vertices)) # all vertices share the same face with the first one, but not cut by the boundary face
        else
            vid = findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x,n)<-EPS,eachcol(vertices)) # all vertices share the same face with the first one, but not cut by the boundary face
        end
        A = hcat(vertices[:,vid],points[pid]...)
        center = vec(mean(A,dims=2))
        id1 = dir%3+1;id2 = (dir+1)%3+1
        phi=@views [atan(x[id2]-center[id2],x[id1]-center[id1]) for x in eachcol(A)]
        id_gauss = sortperm(phi)
        H = 2*(midpoint[dir]-vertices[dir,1])*gaussian_area(@views A[[id1,id2],id_gauss])
        if pos
            return true,8*prod(midpoint-@view(vertices[:,1]))-H,H # gas first
        else
            return true,H,8*prod(midpoint-@view(vertices[:,1]))-H
        end
    else
        local_points = Matrix{Float64}(undef,2,length(points))
        for i in eachindex(points)
            local_points[:,i].=C*points[i]
        end
        center = vec(mean(local_points,dims=2))
        phi=@views [atan(x[2]-center[2],x[1]-center[1]) for x in eachcol(local_points)]
        id = sortperm(phi)
        if length(posid)==4 # isolated vertical face
            centers = @views sum(vertices[:,posid],dims = 2)./4.0
            dir = findfirst(i->abs(centers[i]-midpoint[i])≈0.5*ddu[i],1:3)
            if isnothing(dir)
                H = 0.
            else
                id1 = dir%3+1;id2 = (dir+1)%3+1
                H = pos ? sign(centers[dir]-midpoint[dir])*centers[dir]*ddu[id1]*ddu[id2]/3.0 : -sign(centers[dir]-midpoint[dir])*(2.0*midpoint[dir]-centers[dir])*ddu[id1]*ddu[id2]/3.0
            end
        else
            H = 0.
        end
        if pos
            @views H+=Vertical_Volume_Flux(points[id],midpoint,vertices[:,posid])
            return true,8*prod(midpoint-@view(vertices[:,1]))-H,H # gas first
        else
            negid = findall(x->x<-EPS,dirs)
            @views H+=Vertical_Volume_Flux(points[id],midpoint,vertices[:,negid])
            return true,H,8*prod(midpoint-@view(vertices[:,1]))-H
        end
    end
end