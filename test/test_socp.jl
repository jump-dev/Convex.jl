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
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(sum(square(A*x + b)))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0.42105, TOL)
    @fact evaluate(sum(square(A*x + b))) => roughly(0.42105, TOL)


    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    expr = A * x + b
    p = minimize(sum(expr * expr))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0.42105, TOL)
    @fact evaluate(sum(expr * expr)) => roughly(0.42105, TOL)
  end

  context("inv pos atom") do
    x = Variable(4)
    p = minimize(sum(inv_pos(x)), inv_pos(x) < 2, x > 1, x == 2, 2 == x)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(sum(inv_pos(x))) => roughly(2, TOL)
  end

  context("geo mean atom") do
    x = Variable(2)
    y = Variable(2)
    p = minimize(geo_mean(x, y), x >= 1, y >= 2)
    # not DCP compliant
    @fact vexity(p) => ConcaveVexity()
    p = maximize(geo_mean(x, y), 1 < x, x < 2, y < 2)
    # Just gave it a vector as an objective, not okay
    @fact_throws solve!(p)

    p = maximize(sum(geo_mean(x, y)), 1 < x, x < 2, y < 2)
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(sum(geo_mean(x, y))) => roughly(4, TOL)
  end

  context("sqrt atom") do
    x = Variable()
    p = maximize(sqrt(x), 1 >= x)
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

  context("rational norm atom") do
    A = [-1.175 -1.753  -1.791;
         -0.998 0.446   -0.130;
         1.194  0.978   -1.175];
    B = [0.089  0.617   0.527;
         -0.422 0.596   -1.344;
         -1.650 -0.618  -1.234];
    b = A * ones(3);
    x = Variable(3)
    p = minimize(norm(A * x, 4.5), [B * x == b]);
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(13.218, TOL)
    @fact evaluate(norm(A * x, 4.5)) => roughly(10.9705, TOL)
  end
end
