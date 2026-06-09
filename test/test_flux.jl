using BenchmarkTools, LinearAlgebra
df = rand(10000,2)
vn = rand(10000,2)
@btime @inbounds begin
     for i in axes($df,1)
        for j in axes($df,2)
            $df[i,j]*=$vn[i]
        end
    end
end
@btime begin
    here_micro = [$df[i,j]*$vn[i] for i in axes($df,1),j in axes($df,2)];
end

@btime begin
   for i in axes($df,1) 
        $df[i,:] *= 2.
   end
end

c = zeros(length(df),1)
@btime begin
    @views begin
        for i in axes($df,1)
            $c[i] = dot($df[i,:],$vn[i,:])
        end
    end
end


@btime begin
    @views $df[:,1] .= $df[:,2]
end
@btime begin
    for i in axes($df,1)
        $df[i,1] = $df[i,2]
    end
end

@allocated using KitAMR