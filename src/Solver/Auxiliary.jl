function unpack(t)
    return (getfield(t,i) for i in 1:nfields(t))
end

function MPI_tune(f::Function,args...)
    for i in 1:MPI.Comm_size(MPI.COMM_WORLD)
        if MPI.Comm_rank(MPI.COMM_WORLD)==i-1
            f(args...)
        end
        MPI.Barrier(MPI.COMM_WORLD)
    end
end


function fieldvalues_fn(vs_data,aux_df)
    NDF = typeof(vs_data).parameters[2]
    return [aux_df[:,i] for i in 1:NDF]
end
function fieldvalues_fn(vs_data)
    NDF = typeof(vs_data).parameters[2]
    return [vs_data.df[:,i] for i in 1:NDF]
end
