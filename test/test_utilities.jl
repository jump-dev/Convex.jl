@testset "Utilities" begin

    @testset "clearmemory" begin
        # solve a problem to populate globals
        x = Variable()
        p = minimize(-x, [x <= 0])
        @test vexity(p) == AffineVexity()
        solve!(p, solvers[1])

        @test !isempty(Convex.id_to_variables)
        @test !isempty(Convex.conic_constr_to_constr)

        Convex.clearmemory()

        # check they are cleared
        @test isempty(Convex.id_to_variables)
        @test isempty(Convex.conic_constr_to_constr)
    end

    using SparseArrays
    @testset "Parametrically typed problems" for T = [Float32, Float64, BigFloat]
        x = Variable()
        p = Problem{T}(:minimize, -x, [x <= 0])
        c, A, b, cones, var_to_ranges, vartypes, constraints = Convex.conic_problem(p)
        @test c isa SparseMatrixCSC{T,Int64}
        @test A isa SparseMatrixCSC{T,Int64}
        @test b isa SparseMatrixCSC{T,Int64}
    end

    @testset "ConicObj" for T = [UInt32, UInt64]
        c = ConicObj()
        z = zero(T)
        @test !haskey(c, z)
        c[z] = (1, 1)
        @test c[z] == (1, 1)
        x = T[]
        for (k, v) in c
            push!(x, k)
        end
        @test x == collect(keys(c))
        d = copy(c)
        @test d !== c
    end

    @testset "length and size" begin
        x = Variable(2,3)
        @test length(x) == 6
        @test size(x) == (2, 3)
        @test size(x, 1) == 2
        @test size(x, 2) == 3

        x = Variable(3)
        @test length(x) == 3
        @test size(x) == (3, 1)

        x = Variable()
        @test length(x) == 1
        @test size(x) == (1, 1)
    end

    @testset "lastindex and axes" begin
        x = Variable(2, 3)
        @test axes(x) == (Base.OneTo(2), Base.OneTo(3))
        @test axes(x, 1) == Base.OneTo(2)
        @test lastindex(x) == 6
        @test lastindex(x, 2) == 3
        y = x[:,end]
        @test y isa AbstractExpr
        @test size(y) == (2, 1)
    end

    @testset "Parametric constants" begin
        z = Constant([1.0 0.0im; 0.0 1.0])
        @test z isa Constant{Matrix{Complex{Float64}}}

        # Helper functions
        @test Convex.ispos(1)
        @test Convex.ispos(0)
        @test !Convex.ispos(-1)
        @test Convex.ispos([0,1,0])
        @test !Convex.ispos([0,-1,0])
        @test Convex.isneg(-1)
        @test Convex.isneg(0)
        @test !Convex.isneg(1)
        @test Convex.isneg([0,-1,0])
        @test !Convex.isneg([0,1,0])
        @test Convex._size(3) == (1, 1)
        @test Convex._sign(3) == Positive()
        @test Convex._size([-1,1,1]) == (3, 1)
        @test Convex._sign([-1,1,1]) == NoSign()
        @test Convex._sign([-1,-1,-1]) == Negative()
        @test Convex._size([0 0; 0 0]) == (2, 2)
        @test Convex._sign([0 0; 0 0]) == Positive()
        @test Convex._size(0+1im) == (1, 1)
        @test Convex._sign(0+1im) == ComplexSign()

        @test Convex.imag_conic_form(Constant(1.0)) == [0.0]
        @test Convex.imag_conic_form(Constant([1.0, 2.0])) == [0.0, 0.0]
    end

    @testset "Iteration" begin
        x = Variable(2,3)
        s = sum([xi for xi in x])
        x.value = [1 2 3; 4 5 6]
        # evaluate(s) == [21] (which might be wrong? expected 21)
        # but [21][1] === 21[1] === 21
        # so this should pass even after "fixing" that
        @test evaluate(s)[1] == 21

        x = Variable(4)
        @test [xi.inds for xi in x] == [1:1, 2:2, 3:3, 4:4]

        x = Variable(0)
        @test [xi for xi in x] == []
        @test iterate(x) == nothing
    end
    # returns [21]; not sure why
    # context("iteration") do
    #     x = Variable(2,3)
    #     s = sum([xi for xi in x])
    #     x.value = [1 2 3; 4 5 6]
    #     @fact evaluate(s) --> 21
    # end
end
