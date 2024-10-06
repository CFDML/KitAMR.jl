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
function finalize_p4est!(ps4est::Ptr{p4est_t}, amr::AMR)
    global_data = amr.global_data
    p4est_ghost_destroy(global_data.forest.ghost)
    p4est_mesh_destroy(global_data.forest.mesh)
    pp = PointerWrapper(ps4est)
    p4est_connectivity_destroy(pointer(pp.connectivity))
    p4est_destroy(ps4est)
end
function finalize!(ps4est::Ptr{p4est_t}, amr::AMR)
    finalize_ghost!(amr.ghost.ghost_exchange)
    finalize_IB!(amr.field.IB_buffer)
    finalize_p4est!(ps4est, amr)
end

function finalize_p4est!(ps4est::Ptr{p8est_t}, amr::AMR)
    global_data = amr.global_data
    p8est_ghost_destroy(global_data.forest.ghost)
    p8est_mesh_destroy(global_data.forest.mesh)
    pp = PointerWrapper(ps4est)
    p8est_connectivity_destroy(pointer(pp.connectivity))
    p8est_destroy(ps4est)
end
function finalize_IB!(IB_buffer::Vector{Ptr{Nothing}})
    for i in eachindex(IB_buffer)
        sc_free(-1, IB_buffer[i])
    end
end
function finalize!(ps4est::Ptr{p8est_t}, amr::AMR)
    finalize_ghost!(amr.ghost.ghost_exchange)
    finalize_IB!(amr.field.IB_buffer)
    finalize_p4est!(ps4est, amr)
end
