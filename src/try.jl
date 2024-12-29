using UnPack,BenchmarkTools
struct Test
    a
    b
end
function unpack_f(t::Test)
    return (getfield(t,i) for i in 1:nfields(t))
end
t = Test(1,1)
@btime @unpack a,b=$t
@btime x,y = unpack_f($t)

function test()
    return 1.,2.
end

@btime permutedims(reduce(hcat,[collect(test()) for _ in 1:3]))
@btime begin
    tem = Matrix{Float64}(undef,3,2)
    for i in 1:3
        tem[i,:] .= test()
    end
end

a = collect(1:1000)
b = @views a[200:201]
for i in axes(b)
    @show a[i]
end


