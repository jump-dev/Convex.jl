using Convex, SCS, Test
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
using SparseArrays

TOL = 1e-3

get_solver() = SCS.Optimizer(verbose = 0, eps = 1e-6)

@testset "Some affine problems" begin
    x = Variable(2, 2)
    p = minimize(sum(diag(x, 1)), x >= 1)

    Convex.solve!(p, get_solver())

    @test p.optval ≈ 1 atol = TOL

    x = Variable(10)
    a = rand(10, 1)
    p = maximize(sum(x[2:6]), x <= a)
    @test vexity(p) == AffineVexity()
    solve!(p, get_solver())
    @test p.optval ≈ sum(a[2:6]) atol = TOL
    @test evaluate(sum(x[2:6])) ≈ sum(a[2:6]) atol = TOL


    x = Variable(4, 4)
    p = minimize(sum(diag(x)), x >= 2)
    Convex.solve!(p, get_solver())

    @test p.optval ≈ 8 atol = TOL

end
using LinearAlgebra
@testset "Some sdp tests" begin
    y = Variable((2, 2), :Semidefinite)
    p = minimize(y[1,1])
    Convex.solve!(p, get_solver())
    @test p.optval ≈ 0 atol = TOL


    y = Variable((3, 3))
    p = minimize(sigmamax(y), y[2,1] <= 4, y[2,2] >= 3, sum(y) >= 12)

    E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2)
    s1, s2 = [0.25 -0.25im; 0.25im 0.25], [0.5 0.0; 0.0 0.0]
    p = minimize(real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])
    Convex.solve!(p, get_solver())
    @test p.optval ≈ 1.1464466094719907 rtol = 1e-4
    @test evaluate(tr(E12 * s1 + 2 * E12 * s2 + E21 * s2 + 2 * E21 * s1)) ≈ p.optval rtol = 1e-4
    @test evaluate(E12) ≈ [0.146447 -0.353553im; 0.353553im 0.853553] rtol = 1e-4
    @test evaluate(E21) ≈ [0.853553 0.353553im; -0.353553im 0.146447] rtol = 1e-4
end
