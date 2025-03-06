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
    points = Vector{Float64}[];indices = Int[];i = 1
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
    clp = findall(x->dot(x,n)==0.,vertices)
    if length(clp)+i==1
        return false,0.
    else
        append!(points,vertices)
        append!(indices,CLP)
        sid = sortperm(indices)
        indices = indices[sid]
        points = points[sid]
        aid = findfirst(x->iseven(x)||in(x,CLP[clp]),indices)
        bid = findlast(x->iseven(x)||in(x,CLP[clp]),indices)
        A = Matrix{Float64}(undef,2,bid-aid+1)
        for i in aid:bid
            A[:,i-aid+1] .= points[i]
        end
        l = points[bid]-points[aid]
        if n[1]*l[2]-n[2]*l[1]<0
            return true,gaussian_area(A)/((xmax-xmin)*(ymax-ymin))
        else
            return true,1-gaussian_area(A)/((xmax-xmin)*(ymax-ymin))
        end
    end
end