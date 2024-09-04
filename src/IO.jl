function read_config(filename::T) where {T<:AbstractString}
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @info "reading config from $filename"
        println("--------------------------------------------------------------")
    end
    f = open(filename,"r";lock = false)
    config_dict = Dict{Symbol,Any}()
    for line in eachline(f)
        if length(line) == 0 || line[1] == '#'
            continue
        end
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            println(line)
        end
        var, val = split(line, "=")
        stripped = Symbol(strip(var))
        val = eval(Meta.parse(split(val, "#")[1]))
        config_dict[stripped] = val
    end
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        println("--------------------------------------------------------------")
        println("")
    end
    close(f)
    return config_dict
end