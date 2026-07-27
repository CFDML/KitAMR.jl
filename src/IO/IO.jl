include("Check.jl")
include("Input.jl")
include("Output.jl")
include("Restart.jl")

export read_config,
    save_result, listen_for_save!, check_for_save!, check!, check_for_animsave!
export save_for_restart, restart

function string_to_Cstring(s::String)
    t = codeunits(s)
    p = Ptr{Cstring}(malloc(sizeof(UInt8)*(length(t)+1)))
    ap = Base.unsafe_wrap(Vector{UInt8}, Ptr{UInt8}(p), length(t)+1)
    ap[1:(end-1)] .= t
    ap[end] = 0
    return p
end
