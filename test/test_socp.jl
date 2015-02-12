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
    @fact evaluate(sum_squares(A*x + b))[1] => roughly(0.42105, TOL)
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
    A = [1 2 3; -1 2 3];
    b = A * ones(3);
    x = Variable(3);
    p = minimize(norm(x, 4.5), [A * x == b]);
    @fact vexity(p) => ConvexVexity()
    # Solution is approximately x = [1, .93138, 1.04575]
    solve!(p)
    @fact p.optval => roughly(1.2717, TOL)
    @fact evaluate(norm(x, 4.5)) => roughly(1.2717, TOL)
  end

  context("rational norm dual norm") do
    # Hack to make sure that the SCS Solver runs for 10000 iterations.
    s = get_default_solver()
    if typeof(s).name.name == :SCSSolver
      s = SCSSolver(verbose=0,max_iters=10000);
    end
    v = [0.463339, 0.0216084, -2.07914, 0.99581, 0.889391];
    x = Variable(5);
    q = 1.379;  # q norm constraint that generates many inequalities
    qs = q / (q - 1);  # Conjugate to q
    p = minimize(x' * v);
    p.constraints += (norm(x, q) <= 1);
    @fact vexity(p) => ConvexVexity()
    solve!(p, s)  # Solution is -norm(v, q / (q - 1))
    @fact p.optval => roughly(-2.144087, TOL)
    @fact sum(evaluate(x' * v)) => roughly(-2.144087, TOL)
    @fact evaluate(norm(x, q)) => roughly(1, TOL)
    @fact sum(evaluate(x' * v)) => roughly(-sum(abs(v).^qs)^(1/qs), TOL);
  end

  context("rational norm atom sum") do
    A = [-0.719255  -0.229089;
         -1.33632   -1.37121;
         0.703447  -1.4482];
    b = [-1.82041, -1.67516, -0.866884];
    q = 1.5;
    xvar = Variable(2);
    p = minimize(.5 * sum_squares(xvar) + norm(A * xvar - b, q));
    @fact vexity(p) => ConvexVexity();
    solve!(p)
    # Compute gradient, check it is zero(ish)
    x_opt = xvar.value;
    margins = A * x_opt - b;
    qs = q / (q - 1);  # Conjugate
    denom = sum(abs(margins).^q)^(1/qs);
    g = x_opt + A' * (abs(margins).^(q-1) .* sign(margins)) / denom;
    @fact p.optval => roughly(1.7227, TOL);
    @fact norm(g, 2)^2 => roughly(0, TOL);
  end
end
