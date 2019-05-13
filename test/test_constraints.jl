@testset "Check broadcasted constraints" begin
    R = rand(5,5)
    x = Variable(5)
    y = rand(5)
    # see #218 (https://github.com/JuliaOpt/Convex.jl/issues/218)
    # if this doesn't error, check that it's doing the right thing, i.e.
    # `R*x .<= y` should generate the same constraint as `R*x <= y`.
    @test_throws MethodError R*x .<= y

    p = satisfy(x .== y)
    solve!(p, solvers[1])
    @test x.value ≈ y atol = TOL

    p = satisfy([x == y])
    solve!(p, solvers[1])
    @test x.value ≈ y atol = TOL

end