using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("configure_2D.txt")
ps4est,amr = KitAMR.init(config);

function AMR_4est_corner_iterate(
    forest::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    user_data::Ptr{Nothing},
    corner_iter_fn::Function,
    )
    GC.@preserve forest ghost user_data corner_iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        C_NULL,
        eval(@cfunction($corner_iter_fn, Cvoid, (Ptr{p4est_iter_corner_info}, Ptr{Nothing}))),
    )
end
function AMR_corner_iterate(
    info::Ptr{p4est_iter_corner_info_t},
    data::Ptr{Nothing},
    ::Type{T},
    kernel::Function,
) where {T}
    GC.@preserve info data begin
        ip = PointerWrapper(info)
        GC.@preserve ip kernel(ip, data)
    end
    return nothing
end
function test_corner(info,data)
    AMR_corner_iterate(info, data, KitAMR.P4est_PS_Data, test_corner_kernel)
end
function test_corner_kernel(ip,data)
    d = unsafe_pointer_to_objref(data)
    d[1]+=1
    @show ip.sides.elem_count[]
end

ghost = amr.global_data.forest.ghost
data = [0]
user_data = pointer_from_objref(data)
GC.@preserve data AMR_4est_corner_iterate(ps4est,ghost,user_data,test_corner)
@show data
fp = PointerWrapper(ps4est)
fp.local_num_quadrants[]