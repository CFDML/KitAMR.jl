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
function gaussian_area(A::Matrix) # DIMxN
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
        return false,0.,0.,Float64[],Float64[]
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
            # try B[:,i+1] .= points[index]
            # catch
            #     @show clp bid aid i index length(indices)
            # end
            B[:,i+1] .= points[index]
        end
        l = points[bid]-points[aid]
        if n[1]*l[2]-n[2]*l[1]<0
            solid_weight = gaussian_area(A)
            return true,solid_weight,((xmax-xmin)*(ymax-ymin))-solid_weight,vec(sum(A,dims=2))./size(A,2),vec(sum(B,dims=2))./size(B,2) # solid first
        else
            gas_weight = gaussian_area(A)
            return true,((xmax-xmin)*(ymax-ymin))-gas_weight,gas_weight,vec(sum(B,dims=2))./size(B,2),vec(sum(A,dims=2))./size(A,2)
        end
    end
end
function ray_casting(point::Vector{Float64},vertices::Vector{Vector{Float64}})
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
function find_horizontal_intersection(s_midpoint, f_midpoint, points)
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

function find_vertical_intersection(s_midpoint, f_midpoint, points)
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

function find_intersections(s_midpoint, f_midpoint, closed_curve)
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