function unpack(t)
    return (getfield(t,i) for i in 1:nfields(t))
end