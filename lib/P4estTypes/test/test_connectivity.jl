using MPI
using MPIPreferences
using P4estTypes
using Test

if haskey(ENV, "P4ESTTYPES_TEST_BINARY")
    @test ENV["P4ESTTYPES_TEST_BINARY"] == MPIPreferences.binary
end

MPI.Init()

let
    @test isvalid(Connectivity{4}(:star))
    @test isvalid(Connectivity{8}(:twotrees, 1, 2, 3))
    @test isvalid(Connectivity{4}("data/hole_2d_gmsh.inp"))
    @test isvalid(brick((2, 3)))
    @test isvalid(brick((3, 2, 3)))
    @test brick((3, 4)) == brick((3, 4))
    @test brick((3, 4, 5), (false, false, true)) != brick((3, 4, 5))
    @test isvalid(refine(Connectivity{4}(:star), 2))
    @test isvalid(refine(Connectivity{8}(:twotrees, 1, 2, 3), 2))

    for b in (brick((4, 3)), brick((2, 3, 1)))
        @test sizeof(b) > 0
        P4estTypes.reduce!(b)
        @test isvalid(b)
        P4estTypes.complete!(b)
        @test isvalid(b)
    end
end

MPI.Finalize()
