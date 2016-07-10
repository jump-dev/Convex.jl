using Convex
using FactCheck

facts("Utilities") do

  context("length and size") do
    x = Variable(2,3)
    @fact length(x) --> 6
    @fact size(x) --> (2,3)
    @fact size(x,1) --> 2
    @fact size(x,2) --> 3

    x = Variable(3)
    @fact length(x) --> 3
    @fact size(x) --> (3,1)

    x = Variable()
    @fact length(x) --> 1
    @fact size(x) --> (1,1)
  end

  # returns [21]; not sure why
  # context("iteration") do
  #   x = Variable(2,3)
  #   s = sum([xi for xi in x])
  #   x.value = [1 2 3; 4 5 6]
  #   @fact evaluate(s) --> 21
  # end

end
