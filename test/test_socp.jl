using Convex
using FactCheck

TOL = 1e-3

facts("SOCP Atoms") do

  context("norm 2 atom") do
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(norm_2(A * x + b))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0.64888, TOL)
    @fact evaluate(norm_2(A * x + b)) => roughly(0.64888, TOL)

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    lambda = 1
    p = minimize(norm_2(A * x + b) + lambda * norm_2(x), x >= 1)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(14.9049, TOL)
    @fact evaluate(norm_2(A * x + b) + lambda * norm_2(x)) => roughly(14.9049, TOL)

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    lambda = 1
    p = minimize(norm_2(A * x + b) + lambda * norm_1(x), x >= 1)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(15.4907, TOL)
    @fact evaluate(norm_2(A * x + b) + lambda * norm_1(x)) => roughly(15.4907, TOL)
  end

  context("frobenius norm atom") do
    m = Variable(4, 5)
    c = [m[3, 3] == 4, m >= 1]
    p = minimize(norm(m, :fro), c)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(sqrt(35), TOL)
    @fact evaluate(norm(m, :fro)) => roughly(sqrt(35), TOL)
  end

  context("quad over lin atom") do
    x = Variable(3, 1)
    A = [2 -3 5; -2 9 -3; 5 -8 3]
    b = [-3; 9; 5]
    c = [3 2 4]
    d = -3
    p = minimize(quad_over_lin(A*x + b, c*x + d))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(17.7831, TOL)
    @fact evaluate(quad_over_lin(A*x + b, c*x + d))[1] => roughly(17.7831, TOL)
  end

  context("qol elementwise atom") do
  end

  context("sum squares atom") do
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(sum_squares(A*x + b))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0.42105, TOL)
    @fact evaluate(sum_squares(A*x + b)) => roughly(0.42105, TOL)
  end

  context("square atom") do
  end

  context("inv pos atom") do
  end

  context("geo mean atom") do
  end

  context("sqrt atom") do
  end

  context("quad form atom") do
    x = Variable(3, 1)
    A = [0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
    p = minimize(quad_form(x, A), [x >= 1])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(6.1464, TOL)
    @fact evaluate(quad_form(x, A))[1] => roughly(6.1464, TOL)

    x = Variable(3, 1)
    A = -1.0*[0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
    c = [3 2 4]
    p = maximize(c*x , [quad_form(x, A) >= -1])
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(3.7713, TOL)
    @fact evaluate(quad_form(x, A))[1] => roughly(-1, TOL)
  end

  context("huber atom") do
    x = Variable(3)
    p = minimize(sum(huber(x, 1)), x >= 2)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(9, TOL)
    @fact evaluate(sum(huber(x, 1))) => roughly(9, TOL)
  end

end
