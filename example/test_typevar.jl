
const ts = [Float64,Int64]
struct Test{N}
    x::Float64
end
function test_cfn(p4est::Ptr{p8est_t})
    return nothing
end
function test(t::Test{N}) where{N}
    return @cfunction(test_cfn, Cvoid, (Ptr{ts[N-1]},))
end
test(Test{3}(2.0))