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

function string_to_Cstring(s::String)
    t = codeunits(s)
    p = Ptr{Cstring}(sc_malloc(-1,sizeof(UInt8)*(length(t)+1)))
    ap = Base.unsafe_wrap(Vector{UInt8},Ptr{UInt8}(p),length(t)+1)
    ap[1:end-1] .= t
    ap[end] = 0
    return p
end

#=
example:
function fieldvalues_fn(ps_data)
    return [1/ps_data.prim[end]]
end

KitAMR.write_VTK(ps4est,"testT",["T"],fieldvalues_fn)
=#
function write_VTK(ps4est::Ptr{p8est_t},filename::String,fieldnames::Vector{String},fieldvalues_fn)
    p4est_geo = p8est_geometry_new_connectivity(pointer(PointerWrapper(ps4est).connectivity))
    cont = p8est_vtk_context_new(ps4est,filename)
    p8est_vtk_context_set_geom(cont,p4est_geo)
    p8est_vtk_context_set_continuous(cont,1)
    cont = p8est_vtk_write_header(cont)
    fp = PointerWrapper(ps4est)
    num_quads = fp.local_num_quadrants[]
    ps = Vector{Ptr{Nothing}}(undef,length(fieldnames))
    pas = Vector{Vector{Cdouble}}(undef,length(fieldnames))
    pscs = Vector{Ptr{sc_array_t}}(undef,length(fieldnames))
    for i in eachindex(fieldnames)
        ps[i] = sc_malloc(-1,num_quads*8)
        pas[i] = unsafe_wrap(Vector{Cdouble},Ptr{Cdouble}(ps[i]),num_quads)
        pscs[i] = sc_array_new_data(ps[i],8,num_quads)
    end
    function init_cell_data_kernel(ip,data,dp)
        pas = unsafe_pointer_to_objref(data)
        qid = global_quadid(ip)+1
        ps_data = Base.unsafe_pointer_to_objref(pointer(dp.ps_data))
        fieldvalues = fieldvalues_fn(ps_data)
        for i in eachindex(fieldvalues)
            pas[i][qid] = fieldvalues[i]
        end
    end
    function init_cell_data(info,data)
        AMR_volume_iterate(info, data, P4est_PS_Data, init_cell_data_kernel)
    end
    ppas = pointer_from_objref(pas)
    cfn = @cfunction($init_cell_data, Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})) # explicit cfn to avoid closure error on Windows and MacOS platforms.
    GC.@preserve ps4est ppas pas p8est_iterate(
        ps4est,
        C_NULL,
        ppas,
        cfn,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    pairs = Vector{Union{Ptr{Cstring},Ptr{sc_array_t}}}(undef,2*length(fieldnames))
    p_names = Vector{Ptr{Cstring}}(undef,length(fieldnames))
    for i in eachindex(fieldnames)
        pairs[2*i-1] = p_names[i] = string_to_Cstring(fieldnames[i])
        pairs[2*i] = pscs[i]
    end
    cont = p8est_vtk_write_cell_dataf(cont,1,1,1,0,1,0,pairs...,cont)
    p8est_vtk_write_footer(cont)
    for i in eachindex(fieldnames)
        sc_free(-1,Ptr{Nothing}(p_names[i]))
        sc_free(-1,ps[i])
    end
end

#=
example:
function fieldvalues_fn(vs_data::VS_Data{3})
    return [vs_data.df[:,1],vs_data.level]
end
=#
function write_vs_VTK(vs_data::AbstractVsData{3,1},amr::AMR{3,1},filename::String,fieldnames::Vector{String},fieldvalues_fn)
    global_data = amr.global_data
    xmin,xmax,ymin,ymax,zmin,zmax = global_data.config.quadrature
    Nx,Ny,Nz = global_data.config.vs_trees_num
    AMR_VS_MAXLEVEL = global_data.config.solver.AMR_VS_MAXLEVEL
    dx = (xmax - xmin) / Nx/2^AMR_VS_MAXLEVEL
    dy = (ymax - ymin) / Ny/2^AMR_VS_MAXLEVEL
    dz = (zmax - zmin) / Nz/2^AMR_VS_MAXLEVEL
    D = [dx,dy,dz]
    vertices = Matrix{Float64}(undef,3,8*length(vs_data.level))
    midpoint = vs_data.midpoint
    level = vs_data.level
    dlevel = -(vs_data.level .- AMR_VS_MAXLEVEL)
    cells = Vector{MeshCell}(undef,length(level))
    for i in eachindex(level)
        for j in 1:8
            @. vertices[:,(i-1)*8+j] = midpoint[i,:]+RMT[3][j]*2^dlevel[i]*D/2
        end
        cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL,(1:8).+8*(i-1))
    end
    vtk_grid(filename,vertices,cells;append=false) do vtk
        cell_datas = fieldvalues_fn(vs_data)
        for i in eachindex(fieldnames)
            vtk[fieldnames[i],VTKCellData()] = cell_datas[i]
        end
    end
end