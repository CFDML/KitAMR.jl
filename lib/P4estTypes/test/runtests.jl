using Aqua
using MPI
using MPIPreferences
using P4estTypes
using Pkg
using Test

Aqua.test_all(P4estTypes)

nprocs_str = get(ENV, "P4ESTTYPES_TEST_NPROCS", "")
nprocs = nprocs_str == "" ? clamp(Sys.CPU_THREADS, 2, 4) : parse(Int, nprocs_str)

test_dir = @__DIR__
istest(f) = endswith(f, ".jl") && startswith(f, "test_")
testfiles = sort(filter(istest, readdir(test_dir)))

mktempdir() do tmp_dir
    base_dir = joinpath(@__DIR__, "..")

    # Change to temporary directory so that any files created by the
    # example get cleaned up after execution.
    cd(tmp_dir)
    test_project = Pkg.Types.projectfile_path(test_dir)
    tmp_project = Pkg.Types.projectfile_path(tmp_dir)
    cp(test_project, tmp_project)

    # Copy data files to temporary directory
    test_data_dir = joinpath(test_dir, "data")
    tmp_data_dir = joinpath(tmp_dir, "data")
    mkdir(tmp_data_dir)
    for f in readdir(test_data_dir)
        cp(joinpath(test_data_dir, f), joinpath(tmp_data_dir, f))
    end

    # Setup MPI and P4est preferences
    code = "import Pkg; Pkg.develop(path=raw\"$base_dir\"); Pkg.instantiate(); Pkg.precompile(); include(joinpath(raw\"$test_dir\", \"configure_packages.jl\"))"
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$tmp_project -e "$code"`
    @info "Initializing MPI and P4est with" cmd
    @test success(pipeline(cmd, stderr = stderr, stdout = stdout))

    @testset "$f" for f in testfiles
        cmd = `$(mpiexec()) -n $nprocs $(Base.julia_cmd()) --startup-file=no --project=$tmp_project $(joinpath(test_dir, f))`
        @test success(pipeline(cmd, stderr = stderr, stdout = stdout))
    end
end
