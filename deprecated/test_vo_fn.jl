function AMR_volume_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL,data_type = P4est_PS_Data)
    function iter_fn(info,data)
        GC.@preserve info data begin
            ip = PointerWrapper(info)
            dp = PointerWrapper(data_type, ip.quad.p.user_data[])
            GC.@preserve ip f(ip, data, dp)
        end
    end
    GC.@preserve forest ghost user_data iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        @cfunction($iter_fn,Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
    )
end
function AMR_corner_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL)
    function iter_fn(info,data)
        GC.@preserve info data begin
            ip = PointerWrapper(info)
            GC.@preserve ip f(ip, data)
        end
    end
    GC.@preserve forest ghost user_data iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        C_NULL,
        @cfunction($iter_fn,Cvoid, (Ptr{p4est_iter_corner_info}, Ptr{Nothing})),
    )
end
function AMR_face_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL)
    function iter_fn(info,data)
        GC.@preserve info data begin
            ip = PointerWrapper(info)
            GC.@preserve ip f(ip, data)
        end
    end
    GC.@preserve forest ghost user_data iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        @cfunction($iter_fn,Cvoid, (Ptr{p4est_iter_face_info}, Ptr{Nothing})),
        C_NULL,
    )
end