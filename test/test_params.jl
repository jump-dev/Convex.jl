using Convex
using Base.Test

TOL = 1e-3

@testset "Fixed and freed variables" begin

  @testset "fix and free addition" begin
	x = Variable()
	y = Variable()

	p = minimize(x+y, x>=0, y>=0)
	solve!(p)
	@test isapprox(p.optval, 0, atol=TOL)

	y.value = 4
	fix!(y)
	solve!(p)
	@test isapprox(p.optval, 4, atol=TOL)

	free!(y)
	solve!(p)
	@test isapprox(p.optval, 0, atol=TOL)
  end

  @testset "fix multiplication" begin
	a = [1,2,3,2,1]
	x = Variable(length(a))
	gamma = Variable(Positive())
	fix!(gamma, 0.7)

	p = minimize(norm(x-a) + gamma*norm(x[1:end-1] - x[2:end]))
	solve!(p)
	o1 = p.optval
	# x should be very close to a
	@test isapprox(o1, 0.7 * norm(a[1:end - 1] - a[2:end]), atol=TOL)
	# increase regularization
	fix!(gamma, 1.0)
	solve!(p)
	o2 = p.optval
	# x should be very close to mean(a)
	@test isapprox(o2, norm(a - mean(a)), atol=TOL)

	@test o1 <= o2
  end
end