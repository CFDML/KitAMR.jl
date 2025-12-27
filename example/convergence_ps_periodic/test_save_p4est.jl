using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("./example/periodic/configure_periodic.txt")
ps4est,amr = KitAMR.init(config);
dir = "./example/profile"
pro_dir = pwd()
cd(dir)
p4est_save_ext("p",ps4est,Cint(0),Cint(0))
cd(pro_dir)
@show pwd()
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

function test(;kwargs...)
    if haskey(kwargs,:t)
        @show kwargs[:t]
    end
end