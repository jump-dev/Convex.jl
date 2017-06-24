using Convex
using FactCheck

TOL = 1e-3

facts("Broadcast") do

  x = Variable(5)
  y = [1,2,3,4,5]
  x.value = [1,2,3,4,5]

  context(".*") do
    @fact evaluate(x .* y) --> roughly([1,4,9,16,25], TOL)
    @fact evaluate(y .* x) --> roughly([1,4,9,16,25], TOL)
  end

  context("./") do
    @fact evaluate(x ./ y) --> roughly(ones(Float64, 5), TOL)
    @fact vec(evaluate(y ./ x)) --> roughly(ones(Float64, 5), TOL)
  end

  context(".^") do
    @fact vec(evaluate(x .^ 2)) --> roughly([1.,4.,9.,16.,25.], TOL)
  end

end
