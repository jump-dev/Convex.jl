@testset "Deprecations" begin
    A = Semidefinite(2)
    @test_deprecated lambdamin(A)
    @test_deprecated lambdamax(A)

    @test_deprecated Convex.clearmemory()
end
