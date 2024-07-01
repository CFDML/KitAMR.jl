function finalize_ghost!(ghost_exchange::Ghost_Exchange)
    sc_free(-1, ghost_exchange.ghost_datas)
    sc_free(-1, ghost_exchange.ghost_slopes)
    sc_free(-1, ghost_exchange.ghost_structures)
    for i in eachindex(ghost_exchange.mirror_data_pointers)
        sc_free(-1, ghost_exchange.mirror_data_pointers[i])
        sc_free(-1, ghost_exchange.mirror_slope_pointers[i])
        sc_free(-1, ghost_exchange.mirror_structure_pointers[i])
    end
end
function finalize_p4est!(ps4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    p4est_ghost_destroy(global_data.ghost)
    p4est_mesh_destroy(global_data.mesh)
    pp = PointerWrapper(ps4est)
    p4est_connectivity_destroy(pointer(pp.connectivity))
    p4est_destroy(ps4est)
end
function finalize!(ps4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    finalize_ghost!(DVM_data.ghost_exchange)
    finalize_p4est!(ps4est, DVM_data)
end
