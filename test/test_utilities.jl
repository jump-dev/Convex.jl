using Convex
using Test

@testset "Utilities" begin

  @testset "length and size" begin
    x = Variable(2,3)
    @test length(x) == 6
    @test size(x) == (2, 3)
    @test size(x, 1) == 2
    @test size(x, 2) == 3

    x = Variable(3)
    @test length(x) == 3
    @test size(x) == (3, 1)

    x = Variable()
    @test length(x) == 1
    @test size(x) == (1, 1)
  end

  # returns [21]; not sure why
  # context("iteration") do
  #   x = Variable(2,3)
  #   s = sum([xi for xi in x])
  #   x.value = [1 2 3; 4 5 6]
  #   @fact evaluate(s) --> 21
  # end

end