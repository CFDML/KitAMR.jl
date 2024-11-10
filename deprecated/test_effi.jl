using BenchmarkTools,KitAMR

u = rand(3000);v = rand(3000);midpoint = rand(3000,2)
h = rand(3000);b = rand(3000);df = rand(3000,2)
weight = rand(3000)
@btime KitAMR.micro_to_macro_2D2F($u,$v,$h,$b,$weight)
@btime KitAMR.micro_to_macro_2D2F(@view($midpoint[:,1]),@view($midpoint[:,2]),@view($df[:,1]),@view($df[:,2]),$weight)