using Test

@testset "cut_rect" begin
    rect(xmin,xmax,ymin,ymax) = [[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]]
    as_vertices(vertices) = [collect(Float64,v) for v in vertices]
    run_cut(n,vertices) = KitAMR.cut_rect(collect(Float64,n),as_vertices(vertices))
    rect_area(vertices) = KitAMR.gaussian_area(hcat(as_vertices(vertices)...))

    flag,gas,solid = run_cut([1,1],rect(-1,1,-1,1))
    @test flag
    @test gas ≈ 2.0
    @test solid ≈ 2.0

    flag,gas,solid = run_cut([1,1],rect(0,1,-1,1))
    @test flag
    @test gas ≈ 0.5
    @test solid ≈ 1.5

    flag,gas,solid = run_cut([1,0],rect(-1,1,-1,1))
    @test flag
    @test gas ≈ 2.0
    @test solid ≈ 2.0

    @test run_cut([1,1],rect(0,1,0,1))[1] == false
    @test run_cut([1,1],rect(-1,0,-1,0))[1] == false

    for n in ([1.0,0.3],[-0.7,1.2],[1.0,1.0+1e-13],[1e-7,1.0])
        for vertices in (rect(-1,1,-0.5,0.75),rect(-0.2,1.7,-1.3,0.4),rect(-3,-1,1,2))
            flag,gas,solid = run_cut(n,vertices)
            if flag
                @test gas > 0.0
                @test solid > 0.0
                @test gas+solid ≈ rect_area(vertices)
            end
        end
    end
end

@testset "cut_cube" begin
    function cube_vertices(midpoint,ddu)
        vertices = Matrix{Float64}(undef,3,8)
        for j in 1:8
            vertices[:,j] .= 0.5 .* KitAMR.ANTIVT[3][j] .* ddu .+ midpoint
        end
        return vertices
    end

    n = [1.0,1.0,1.0]
    n ./= sqrt(sum(abs2,n))
    midpoint = [0.0,0.0,0.0]
    ddu = [1.0,1.0,1.0]
    flag,gas,solid = KitAMR.cut_cube(n,KitAMR.cut_cube_rotate(n),midpoint,ddu,cube_vertices(midpoint,ddu))
    @test flag
    @test gas ≈ 0.5
    @test solid ≈ 0.5
end
