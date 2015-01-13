using Convex
using FactCheck

TOL = 1e-3

facts("LP Atoms") do

  context("abs atom") do
    x = Variable()
    p = minimize(abs(x), x<=-1)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)
    @fact evaluate(abs(x)) => roughly(1, TOL)

    x = Variable(2,2)
    p = minimize(sum(abs(x)), x[2,2]>=1, x[1,1]>=1, x>=0)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(sum(abs(x))) => roughly(2, TOL)
  end

  context("maximum atom") do
    x = Variable(10)
    a = rand(10, 1)
    p = minimize(maximum(x), x >= a)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(maximum(a), TOL)
    @fact evaluate(maximum(x)) => roughly(maximum(a), TOL)
  end

  context("minimum atom") do
    x = Variable(1)
    a = rand(10, 10)
    p = maximize(minimum(x), x <= a)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(minimum(a), TOL)
    @fact evaluate(minimum(x)) => roughly(minimum(a), TOL)

    x = Variable(4, 4)
    y = Variable(4, 6)
    z = Variable(1)
    c = ones(4, 1)
    d = 2 * ones(6, 1)
    constraints = [[x y] <= 2, z <= 0, z <= x, 2z >= -1]
    objective = sum(x + z) + minimum(y) + c' * y * d
    p = maximize(objective, constraints)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(130, TOL)
    @fact evaluate(objective)[1] => roughly(130, TOL)
  end

  context("max atom") do
    x = Variable(10, 10)
    y = Variable(10, 10)
    a = rand(10, 10)
    b = rand(10, 10)
    p = minimize(maximum(max(x, y)), [x >= a, y >= b])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    max_a = maximum(a)
    max_b = maximum(b)
    @fact p.optval => roughly(max(max_a, max_b), TOL)
    @fact evaluate(maximum(max(x, y))) => roughly(max(max_a, max_b), TOL)
  end

  context("min atom") do
    x = Variable(10, 10)
    y = Variable(10, 10)
    a = rand(10, 10)
    b = rand(10, 10)
    p = maximize(minimum(min(x, y)), [x <= a, y <= b])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    min_a = minimum(a)
    min_b = minimum(b)
    @fact p.optval => roughly(min(min_a, min_b), TOL)
    @fact evaluate(minimum(min(x, y))) => roughly(min(min_a, min_b), TOL)
  end

  context("pos atom") do
    x = Variable(3)
    a = [-2; 1; 2]
    p = minimize(sum(pos(x)), [x >= a, x <= 2])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)
    @fact evaluate(sum(pos(x))) => roughly(3, TOL)
  end

  context("neg atom") do
  end

  context("hinge loss atom") do
  end

  context("norm inf atom") do
    x = Variable(3)
    p = minimize(norm_inf(x), [-2 <= x, x <= 1])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(norm_inf(x)) => roughly(0, TOL)
  end

  context("norm 1 atom") do
    x = Variable(3)
    p = minimize(norm_1(x), [-2 <= x, x <= 1])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(norm_1(x)) => roughly(0, TOL)
  end

end
