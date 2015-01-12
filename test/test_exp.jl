using Convex
using FactCheck

TOL = 1e-3

facts("Exp Atoms") do

  context("exp atom") do
    y = Variable()
    p = minimize(exp(y), y>=0)
    solve!(p)
    @fact p.optval => roughly(1, TOL)

    y = Variable()
    p = minimize(exp(y), y>=1)
    solve!(p)
    @fact p.optval => roughly(exp(1), TOL)

    y = Variable(5)
    p = minimize(sum(exp(y)), y>=0)
    solve!(p)
    @fact p.optval => roughly(5, TOL)

    y = Variable(5)
    p = minimize(sum(exp(y)), y>=0)
    solve!(p)
    @fact p.optval => roughly(5, TOL)
  end

  context("log atom") do
    y = Variable()
    p = maximize(log(y), y<=1)
    solve!(p)
    @fact p.optval => roughly(0, TOL)

    y = Variable()
    p = maximize(log(y), y<=2)
    solve!(p)
    @fact p.optval => roughly(log(2), TOL)

    y = Variable()
    p = maximize(log(y), [y<=2, exp(y)<=10])
    solve!(p)
    @fact p.optval => roughly(log(2), TOL)
  end

  context("log sum exp atom") do
    y = Variable(5);
    p = minimize(logsumexp(y), y>=1);
    solve!(p)
    @fact p.optval => roughly(log(exp(1)*5), TOL)
  end

  context("entropy atom") do
    y = Variable(5, Positive());
    p = maximize(entropy(y), sum(y)<=1);
    solve!(p)
    @fact p.optval => roughly(-log(1/5), TOL)
  end

end
