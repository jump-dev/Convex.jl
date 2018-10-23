using Convex
using Test

TOL = 1e-3

@testset "SOCP Atoms" begin

  @testset "norm 2 atom" begin
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(norm2(A * x + b))
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0.64888, atol=TOL)
    @test isapprox(evaluate(norm2(A * x + b)), 0.64888, atol=TOL)

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    lambda = 1
    p = minimize(norm2(A * x + b) + lambda * norm2(x), x >= 1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 14.9049, atol=TOL)
    @test isapprox(evaluate(norm2(A * x + b) + lambda * norm2(x)), 14.9049, atol=TOL)

    x = Variable(2)

    p = minimize(norm2([x[1] + 2x[2] + 2; 2x[1] + x[2] + 3; 3x[1]+4x[2] + 4]) + lambda * norm2(x), x >= 1)
    @test vexity(p) == ConvexVexity()

    solve!(p)
    @test isapprox(p.optval, 14.9049, atol=TOL)
    @test isapprox(evaluate(norm2(A * x + b) + lambda * norm2(x)), 14.9049, atol=TOL)

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    lambda = 1
    p = minimize(norm2(A * x + b) + lambda * norm_1(x), x >= 1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 15.4907, atol=TOL)
    @test isapprox(evaluate(norm2(A * x + b) + lambda * norm_1(x)), 15.4907, atol=TOL)
  end

# TODO: fix ambiguity
# @testset "frobenius norm atom" begin
#   m = Variable(4, 5)
#   c = [m[3, 3] == 4, m >= 1]
#   p = minimize(vecnorm(m, 2), c)
#   @test vexity(p) == ConvexVexity()
#   solve!(p)
#   @test isapprox(p.optval, sqrt(35), atol=TOL)
#   @test isapprox(evaluate(vecnorm(m, 2)), sqrt(35), atol=TOL)
# end

  @testset "quad over lin atom" begin
    x = Variable(3, 1)
    A = [2 -3 5; -2 9 -3; 5 -8 3]
    b = [-3; 9; 5]
    c = [3 2 4]
    d = -3
    p = minimize(quadoverlin(A*x + b, c*x + d))
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 17.7831, atol=TOL)
# TODO: fix iterate(::Convex.MultiplyAtom)
#   @test isapprox((evaluate(quadoverlin(A * x + b, c * x + d)))[1], 17.7831, atol=TOL)
  end

  @testset "sum squares atom" begin
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(sumsquares(A*x + b))
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0.42105, atol=TOL)
    @test isapprox((evaluate(sumsquares(A * x + b)))[1], 0.42105, atol=TOL)
  end

  @testset "square atom" begin
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(sum(square(A*x + b)))
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0.42105, atol=TOL)
    @test isapprox(evaluate(sum(square(A * x + b))), 0.42105, atol=TOL)

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    expr = A * x + b
    p = minimize(sum(dot(^)(expr,2))) # elementwise ^
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0.42105, atol=TOL)
    @test isapprox(evaluate(sum(broadcast(^, expr, 2))), 0.42105, atol=TOL)

    p = minimize(sum(dot(*)(expr, expr))) # elementwise *
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0.42105, atol=TOL)
    @test isapprox(evaluate(sum((dot(*))(expr, expr))), 0.42105, atol=TOL)
  end

  @testset "inv pos atom" begin
    x = Variable(4)
    p = minimize(sum(invpos(x)), invpos(x) < 2, x > 1, x == 2, 2 == x)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 2, atol=TOL)
    @test isapprox(evaluate(sum(invpos(x))), 2, atol=TOL)

    x = Variable(3)
    p = minimize(sum(dot(/)([3,6,9], x)), x<=3)
    solve!(p)
    @test isapprox(x.value, 3 * ones(3, 1), atol=TOL)
    @test isapprox(p.optval, 6, atol=TOL)
    @test isapprox(evaluate(sum((dot(/))([3, 6, 9], x))), 6, atol=TOL)

    x = Variable()
    p = minimize(sum([3,6,9]/x), x<=3)
    solve!(p)
    @test isapprox(x.value, 3, atol=TOL)
    @test isapprox(p.optval, 6, atol=TOL)
    @test isapprox(evaluate(sum([3, 6, 9] / x)), 6, atol=TOL)
  end

  @testset "geo mean atom" begin
    x = Variable(2)
    y = Variable(2)
    p = minimize(geomean(x, y), x >= 1, y >= 2)
    # not DCP compliant
    @test vexity(p) == ConcaveVexity()
    p = maximize(geomean(x, y), 1 < x, x < 2, y < 2)
    # Just gave it a vector as an objective, not okay
    @test_throws Exception solve!(p)

    p = maximize(sum(geomean(x, y)), 1 < x, x < 2, y < 2)
    solve!(p)
    @test isapprox(p.optval, 4, atol=TOL)
    @test isapprox(evaluate(sum(geomean(x, y))), 4, atol=TOL)
  end

  @testset "sqrt atom" begin
    x = Variable()
    p = maximize(sqrt(x), 1 >= x)
  end

  @testset "quad form atom" begin
    x = Variable(3, 1)
    A = [0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
    p = minimize(quadform(x, A), [x >= 1])
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 6.1464, atol=TOL)
    @test isapprox((evaluate(quadform(x, A)))[1], 6.1464, atol=TOL)

    x = Variable(3, 1)
    A = -1.0*[0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
    c = [3 2 4]
    p = maximize(c*x , [quadform(x, A) >= -1])
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 3.7713, atol=TOL)
    @test isapprox((evaluate(quadform(x, A)))[1], -1, atol=TOL)
  end

  @testset "huber atom" begin
    x = Variable(3)
    p = minimize(sum(huber(x, 1)), x >= 2)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 9, atol=TOL)
    @test isapprox(evaluate(sum(huber(x, 1))), 9, atol=TOL)
  end

# @testset "rational norm atom" begin
#   A = [1 2 3; -1 2 3]
#   b = A * ones(3)
#   x = Variable(3)
# TODO: fix norm ambiguity
#   p = minimize(norm(x, 4.5), [A * x == b])
#   @test vexity(p) == ConvexVexity()
    # Solution is approximately x = [1, .93138, 1.04575]
#   solve!(p)
#   @test isapprox(p.optval, 1.2717, atol=TOL)
#   @test isapprox(evaluate(norm(x, 4.5)), 1.2717, atol=TOL)
# end

# @testset "rational norm dual norm" begin
# TODO: fix norm ambiguity
#   v = [0.463339, 0.0216084, -2.07914, 0.99581, 0.889391]
#   x = Variable(5)
#   q = 1.379;  # q norm constraint that generates many inequalities
#   qs = q / (q - 1);  # Conjugate to q
#   p = minimize(x' * v)
#   p.constraints += (norm(x, q) <= 1)
#   @test vexity(p) == ConvexVexity()
#   solve!(p)  # Solution is -norm(v, q / (q - 1))
#   @test isapprox(p.optval, -2.144087, atol=TOL)
#   @test isapprox(sum(evaluate(x' * v)), -2.144087, atol=TOL)
#   @test isapprox(evaluate(norm(x, q)), 1, atol=TOL)
#   @test isapprox(sum(evaluate(x' * v)), -(sum(abs.(v) .^ qs) ^ (1 / qs)), atol=TOL)
# end

# @testset "rational norm atom sum" begin
# TODO: fix norm ambiguity
#   A = [-0.719255  -0.229089
#        -1.33632   -1.37121
#        0.703447  -1.4482]
#   b = [-1.82041, -1.67516, -0.866884]
#   q = 1.5
#   xvar = Variable(2)
#   p = minimize(.5 * sumsquares(xvar) + norm(A * xvar - b, q))
#   @test vexity(p) == ConvexVexity()
#   solve!(p)
    # Compute gradient, check it is zero(ish)
#   x_opt = xvar.value
#   margins = A * x_opt - b
#   qs = q / (q - 1);  # Conjugate
#   denom = sum(abs.(margins).^q)^(1/qs)
#   g = x_opt + A' * (abs.(margins).^(q-1) .* sign.(margins)) / denom
#   @test isapprox(p.optval, 1.7227, atol=TOL)
#   @test isapprox(norm(g, 2) ^ 2, 0, atol=TOL)
# end

# @testset "norm consistent with Base" begin
# TODO: fix norm ambiguity
#   A = randn(4, 4)
#   x = Variable(4, 4)
#   x.value = A
#   @test isapprox(evaluate(norm(x)), norm(A), atol=TOL)
#   @test isapprox(evaluate(norm(x, 1)), norm(A, 1), atol=TOL)
#   @test isapprox(evaluate(norm(x, 2)), norm(A, 2), atol=TOL)
#   @test isapprox(evaluate(norm(x, Inf)), norm(A, Inf), atol=TOL)
#   @test isapprox(evaluate(vecnorm(x, 1)), norm(vec(A), 1), atol=TOL)
#   @test isapprox(evaluate(vecnorm(x, 2)), norm(vec(A), 2), atol=TOL)
#   @test isapprox(evaluate(vecnorm(x, 7)), norm(vec(A), 7), atol=TOL)
#   @test isapprox(evaluate(vecnorm(x, Inf)), norm(vec(A), Inf), atol=TOL)
# end


end
