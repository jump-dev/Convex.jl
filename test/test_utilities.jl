@testset "Utilities" begin

    @testset "ConicObj" begin
        c = ConicObj()
        z = UInt64(0)
        @test !haskey(c, z)
        c[z] = (1, 1)
        @test c[z] == (1, 1)
        x = UInt64[]
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

    # returns [21]; not sure why
    # context("iteration") do
    #     x = Variable(2,3)
    #     s = sum([xi for xi in x])
    #     x.value = [1 2 3; 4 5 6]
    #     @fact evaluate(s) --> 21
    # end
end
