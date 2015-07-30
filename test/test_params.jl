using Convex
using FactCheck

TOL = 1e-3

facts("Fixed and freed variables") do

  context("fix and free addition") do
	x = Variable()
	y = Variable()

	p = minimize(x+y, x>=0, y>=0)
	solve!(p)
	@fact p.optval => roughly(0, TOL)

	y.value = 4
	fix!(y)
	solve!(p)
	@fact p.optval => roughly(4, TOL)

	free!(y)
	solve!(p)
	@fact p.optval => roughly(0, TOL)	
  end

  context("fix multiplication") do	
	a = [1,2,3,2,1]
	x = Variable(length(a))
	gamma = Variable(Positive())
	fix!(gamma, 1)

	p = minimize(norm(x-a) + gamma*norm(x[1:end-1] - x[2:end]))
	solve!(p)
	o1 = p.optval

	# increase regularization
	fix!(gamma, 3)
	solve!(p)
	o2 = p.optval

	@fact o1 < o2 => true
  end  
end
