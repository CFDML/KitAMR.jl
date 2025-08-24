function pause(ps4est,amr;path::String = "pause_data",kwargs...)
    save_pause_data(ps4est,amr,path)
    save_pause_p4est(ps4est;path)
end
function save_pause_data(ps4est,amr::AMR{DIM,NDF},path::String) where{DIM,NDF}
    fp = PointerWrapper(ps4est)
    N = fp.local_num_quadrants[];rank = MPI.Comm_rank(MPI.COMM_WORLD)
    vs_nums = Vector{Int}(undef,N);bound_encs = Vector{Int}(undef,N)
    index = 1
    for tree in amr.field.trees.data
        for ps_data in tree
            if isa(ps_data,InsideSolidData)
                vs_nums[index] = 0
                bound_encs[index] = 0
            else
                vs_nums[index] = ps_data.vs_data.vs_num
                bound_encs[index] = ps_data.bound_enc
            end
            index += 1
        end
    end
    Nv = sum(vs_nums)
    levels = Vector{Int8}(undef,Nv)
    midpoints = Matrix{Float64}(undef,Nv,DIM)
    dfs = Matrix{FLoat64}(undef,Nv,NDF)
    index = 0
    for tree in amr.field.trees.data
        for ps_data in tree
            if !isa(ps_data,InsideSolidData)
                vs_data = ps_data.vs_data
                vs_num = vs_data.vs_num
                range = index+1:index+vs_num
                levels[range].=vs_data.level
                midpoints[range,:].=vs_data.midpoint
                dfs[range,:].=vs_data.df
            end
            index+=1
        end
    end
    jldopen(path*"pause_data"*string(rank)*".jld2","w";compress=true) do f
        f["vs_nums"]=vs_nums
        f["bound_encs"]=bound_encs
        f["levels"]=levels
        f["midpoints"]=midpoints
        f["dfs"]=dfs
    end
    if rank==0
        jldopen(path*"pause_data_set.jld2","w") do f
            f["gfq"]=unsafe_wrap(Vector{Int},pointer(fp.global_first_quadrant),MPI.Comm_size(MPI.COMM_WORLD)+1)
            f["config"]=amr.global_data.config
            f["status"]=StatusForSave(amr.global_data.status)
        end
    end
end
function save_pause_p4est(p4est::Ptr{p4est_t};path::String)
    pro_path = pwd()
    cd(path)
    GC.@preserve p4est p4est_save_ext("p",p4est,Cint(0),Cint(0))
    cd(pro_path)
    return nothing
end
function restart(dirname::String)
    set = load(dirname*"paus_data_set.jld2")
    global_data = recons_global_data(set[:config],set[:status])
    ps4est = pxest_load(dirname,global_data)
    pp = PointerWrapper(ps4est)
    gfq = Base.unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    gfq_ = set[:gfq]
    first_rank = findfirst(x->gfq[rank+1]<=x,gfq_)-1
    last_rank = findlast(x->gfq[rank+2]>x,gfq_)-1
    datas = [load(dirname*"pause_data"*string(x)*".jld2") for x in first_rank:last_rank]
    first_ps_offset = gfq[rank+1]-gfq_[first_rank+1]
    last_ps_offset = gfq[rank+2]-gfq_[last_rank+1]
    vs_nums = datas[1][:vs_nums][first_ps_offset+1:end]
    first_vs_offset = sum(vs_nums)
    bound_encs = datas[1][:bound_encs][first_ps_offset+1:end]
    if length(datas)>1
        for i in 2:length(datas)-1
            append!(vs_nums,datas[i][:vs_nums])
            append!(bound_encs,datas[i][:bound_encs])
        end
        last_vs_nums = datas[end][:vs_nums][1:last_ps_offset]
        last_vs_offset = sum(last_vs_nums)
        append!(vs_nums,last_vs_nums)
        append!(bound_encs,datas[end][:bound_encs][1:last_ps_offset])
    end
    levels = datas[1][:levels][first_vs_offset+1:end]
    midpoints = datas[1][:midpoints][first_vs_offset+1:end,:]
    
end
function recons_global_data(config::Configure{DIM,NDF},stat::StatusForSave) where{DIM,NDF}
    global_data = Global_Data(config)
    global_data.status = Status(stat)
    return global_data
end
function pxest_load(dirname::String,global_data::Global_Data{DIM,NDF};kwargs...) where{DIM,NDF}
    pro_path = pwd()
    cd(dirname)
    if DIM==2
        cnn = Ptr{Ptr{p4est_connectivity_t}}(Libc.malloc(sizeof(Ptr{Ptr{p4est_connectivity_t}})))
        ps4est = GC.@preserve cnn p4est_load_ext("p",MPI.COMM_WORLD,Cint(0),Cint(0),Cint(1),Cint(0),C_NULL,cnn)
    else
        cnn = Ptr{Ptr{p8est_connectivity_t}}(Libc.malloc(sizeof(Ptr{Ptr{p8est_connectivity_t}})))
        ps4est = GC.@preserve cnn p8est_load_ext("p",MPI.COMM_WORLD,Cint(0),Cint(0),Cint(1),Cint(0),C_NULL,cnn)
    end
    cd(pro_path)
    return ps4est
end