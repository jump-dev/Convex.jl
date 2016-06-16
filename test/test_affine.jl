using Convex
using FactCheck

TOL = 1e-3

facts("Affine Atoms") do

  context("negate atom") do
    x = Variable()
    p = minimize(-x, [x <= 0])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(-x) => roughly(0, TOL)
  end

  context("multiply atom") do
    x = Variable(1)
    p = minimize(2.0 * x, [x >= 2, x <= 4])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(2.0 * x)[1] => roughly(4, TOL)

    x = Variable(2)
    A = 1.5 * eye(2)
    p = minimize([2 2] * x, [A * x >= [1.1; 1.1]])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(2.93333, TOL)
    @fact evaluate([2 2] * x)[1] => roughly(2.93333, TOL)
    @fact vec(evaluate(A * x)) => roughly([1.1; 1.1], TOL)

    y = Variable(1)
    x = Variable(3)
    z = [1.0, 2.0, 3.0] * y
    k = -y * [1.0, 2.0, 3.0]
    c = [y <= 3.0, y >= 0.0, x >= ones(3), k <= x, x <= z]
    o = 3 * y
    p = Problem(:minimize, o, c)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)

    p = Problem(:minimize, o, c...)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)
  end

  context("dot atom") do
    x = Variable(2)
    p = minimize(dot([2.0; 2.0], x), x >= [1.1; 1.1])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(4.4, TOL)
    @fact evaluate(dot([2.0; 2.0], x))[1] => roughly(4.4, TOL)
  end

  context("vecdot atom") do
    x = Variable(2,2)
    p = minimize(vecdot(2*ones(2,2), x), x >= 1.1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(8.8, TOL)
    @fact evaluate(vecdot(2*ones(2,2), x))[1] => roughly(8.8, TOL)
  end

  context("add atom") do
    x = Variable(1)
    y = Variable(1)
    p = minimize(x + y, [x >= 3, y >= 2])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(5, TOL)
    @fact evaluate(x + y) => roughly(5, TOL)

    x = Variable(1)
    p = minimize(x, [eye(2) + x >= eye(2)])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(eye(2) + x) => roughly(eye(2), TOL)

    y = Variable()
    p = minimize(y - 5, y >= -1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(-6, TOL)
    @fact evaluate(y - 5) => roughly(-6, TOL)
  end

  context("transpose atom") do
    x = Variable(2)
    c = ones(2, 1)
    p = minimize(x' * c, x >= 1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(x' * c)[1] => roughly(2, TOL)

    X = Variable(2, 2)
    c = ones(2, 1)
    p = minimize(c' * X' * c, [X >= ones(2, 2)])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(c' * X' * c)[1] => roughly(4, TOL)

    rows = 2
    cols = 3
    r = rand(rows, cols)
    r_2 = rand(cols, rows)
    x = Variable(rows, cols)
    c = ones(1, cols)
    d = ones(rows, 1)
    p = minimize(c * x' * d + d' * x * c' + (c * x''''' * d)',
                [x' >= r_2, x >= r, x''' >= r_2, x'' >= r])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    s = sum(max(r, r_2')) * 3
    @fact p.optval => roughly(s, TOL)
    @fact evaluate(c * x' * d + d' * x * c' + (c * x''''' * d)')[1] => roughly(s, TOL)
  end

  context("index atom") do
    x = Variable(2)
    p = minimize(x[1] + x[2], [x >= 1])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(x[1] + x[2])[1] => roughly(2, TOL)

    x = Variable(3)
    I = [true true false]
    p = minimize(sum(x[I]), [x >= 1])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(sum(x[I]))[1] => roughly(2, TOL)

    rows = 6
    cols = 8
    n = 2
    X = Variable(rows, cols)
    A = randn(rows, cols)
    c = rand(1, n)
    p = minimize(c * X[1:n, 5:5+n-1]' * c', X >= A)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    s = c * A[1:n, 5:5+n-1]' * c'
    @fact p.optval => roughly(s[1], TOL)
    @fact evaluate(c * X[1:n, 5:5+n-1]' * c') => roughly(s, TOL)
  end

  context("sum atom") do
    x = Variable(2,2)
    p = minimize(sum(x), x>=1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(sum(x)) => roughly(4, TOL)

    x = Variable(2,2)
    p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)
    @fact evaluate(sum(x) - 2*x[1,1])[1] => roughly(1, TOL)

    x = Variable(10)
    a = rand(10, 1)
    p = maximize(sum(x[2:6]), x <= a)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(sum(a[2:6]), TOL)
    @fact evaluate(sum(x[2:6])) => roughly(sum(a[2:6]), TOL)
  end

  context("diag atom") do
    x = Variable(2,2)
    p = minimize(sum(diag(x,1)), x >= 1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)
    @fact evaluate(sum(diag(x,1))) => roughly(1, TOL)

    x = Variable(4, 4)
    p = minimize(sum(diag(x)), x >= 2)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(8, TOL)
    @fact evaluate(sum(diag(x))) => roughly(8, TOL)
  end

  context("trace atom") do
    x = Variable(2,2)
    p = minimize(trace(x), x >= 1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(trace(x)) => roughly(2, TOL)
  end

  context("dot multiply atom") do
    x = Variable(3)
    p = maximize(sum(x.*[1,2,3]), x<=1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(6, TOL)
    @fact evaluate(sum(x.*[1,2,3])) => roughly(6, TOL)

    x = Variable(3, 3)
    p = maximize(sum(x.*eye(3)), x<=1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)
    @fact evaluate(sum(x.*eye(3))) => roughly(3, TOL)

    x = Variable(5, 5)
    p = minimize(x[1, 1], 3 .* x >= 3)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)
    @fact evaluate(x[1, 1])[1] => roughly(1, TOL)

    x = Variable(1, 3, Positive())
    p = maximize(sum(x./[1 2 3]), x<=1)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(11/6, TOL)
    @fact evaluate(sum(x./[1 2 3])) => roughly(11/6, TOL)
  end

  context("reshape atom") do
    A = rand(2, 3)
    X = Variable(3, 2)
    c = rand()
    p = minimize(sum(reshape(X, 2, 3) + A), X >= c)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(sum(A + c), TOL)
    @fact evaluate(sum(reshape(X, 2, 3) + A)) => roughly(sum(A + c), TOL)

    b = rand(6)
    p = minimize(sum(vec(X) + b), X >= c)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(sum(b + c), TOL)
    @fact evaluate(sum(vec(X) + b)) => roughly(sum(b + c), TOL)

    x = Variable(4, 4)
    c = ones(16, 1)
    reshaped = reshape(x, 16, 1)
    a = collect(1:16)
    p = minimize(c' * reshaped, reshaped >= a)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    # TODO: why is accuracy lower here?
    @fact p.optval => roughly(136, 10*TOL)
    @fact evaluate(c' * reshaped)[1] => roughly(136, 10*TOL)
  end

  context("hcat atom") do
    x = Variable(4, 4)
    y = Variable(4, 6)
    p = maximize(sum(x) + sum([y 4*ones(4)]), [x y 2*ones(4, 2)] <= 2)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(96, TOL)
    @fact evaluate(sum(x) + sum([y 4*ones(4)])) => roughly(96, TOL)
    @fact evaluate([x y 2*ones(4,2)]) => roughly(2*ones(4, 12), TOL)
  end

  context("vcat atom") do
    x = Variable(4, 4)
    y = Variable(4, 6)
    p = maximize(sum(x) + sum([y 4*eye(4); x -ones(4, 6)]), [x, y'] <= 2)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    # TODO: why is accuracy lower here?
    @fact p.optval => roughly(104, 10*TOL)
    @fact evaluate(sum(x) + sum([y 4*eye(4); x -ones(4, 6)])) => roughly(104, 10*TOL)
    @fact evaluate([x, y']) => roughly(2*ones(10, 4), TOL)
  end

  context("diagm atom") do
    x = Variable(2, 2)
    @fact_throws diagm(x)

    x = Variable(4)
    p = minimize(sum(diagm(x)), x == [1; 2; 3; 4])
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(10, TOL)
    @fact all(abs(evaluate(diagm(x)) - diagm([1; 2; 3; 4])) .<= TOL) => true

    x = Variable(3)
    c = [1; 2; 3]
    p = minimize(c' * diagm(x) * c, x >= 1, sum(x) == 10)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(21, TOL)

    x = Variable(3)
    p = minimize(sum(x), x >= 1, diagm(x)[1, 2] == 1)
    @fact solve!(p) => nothing
    @fact (p.status != :Optimal) => true
  end

  context("conv atom") do
    x = Variable(3)
    h = [1, -1]
    p = minimize(sum(conv(h, x)) + sum(x), x >= 1, x <= 2)
    @fact vexity(p) => AffineVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)
    @fact evaluate(sum(conv(h, x))) => roughly(0, TOL)
  end

  context("satisfy problems") do
    x = Variable()
    p = satisfy(x >= 0)
    add_constraints!(p, x >= 1)
    add_constraints!(p, [x >= -1, x <= 4])
    solve!(p)
    @fact p.status => :Optimal

    p = satisfy([x >= 0, x >= 1, x <= 2])
    solve!(p)
    @fact p.status => :Optimal

    p = maximize(1, [x >= 1, x <= 2])
    solve!(p)
    @fact p.status => :Optimal

    constr = x >= 0
    constr += x >= 1
    constr += x <= 10
    constr2 = x >= 0
    constr2 += [x >= 2, x <= 3] + constr
    p = satisfy(constr)
    solve!(p)
    @fact p.status => :Optimal
  end

  context("dual") do
    x = Variable()
    p = minimize(x, x >= 0)
    solve!(p)
    if p.solution.has_dual
        println("Solution object has dual value, checking for dual correctness.")
        @fact p.constraints[1].dual => roughly(1, TOL)
    end

    x = Variable()
    p = maximize(x, x <= 0)
    solve!(p)
    if p.solution.has_dual
        println("Solution object has dual value, checking for dual correctness.")
        @fact p.constraints[1].dual => roughly(1, TOL)
    end

    x = Variable()
    p = minimize(x, x >= 0, x == 2)
    solve!(p)
    if p.solution.has_dual
        println("Solution object has dual value, checking for dual correctness.")
        @fact p.constraints[1].dual => roughly(0, TOL)
        @fact abs(p.constraints[2].dual) => roughly(1, TOL)
    end

    x = Variable(2)
    A = 1.5 * eye(2)
    p = minimize(dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
    solve!(p)
    if p.solution.has_dual
        println("Solution object has dual value, checking for dual correctness.")
        dual = [4/3; 4/3]
        @fact all(abs(p.constraints[1].dual - dual) .<= TOL) => true
    end
  end
end
