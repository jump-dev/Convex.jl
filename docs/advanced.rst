=====================================
Advanced Features
=====================================

Dual Variables
******************

Convex.jl also returns the optimal dual variables for a problem. These are stored in the :code:`dual` field associated with each constraint.
::

	using Convex

	x = Variable()
	constraint = x >= 0
	p = minimize(x, constraint)
	solve!(p)

	# Get the dual value for the constraint
	p.constraints[1].dual
	# or
	constraint.dual

Warmstarting
******************

If you're solving the same problem many times with different values
of a parameter, Convex.jl can initialize many solvers with the solution
to the previous problem, which sometimes speeds up the solution time.
This is called a **warm start**. 

To use this feature,
pass the optional argument `warmstart=true` to the `solve!` method.
::

	# initialize data
	n = 1000
	y = rand(n)
	x = Variable(n)

	# first solve
	lambda = 100
	problem = minimize(sumsquares(y - x) + lambda * sumsquares(x - 10))
	@time solve!(problem)

	# now warmstart
	# if the solver takes advantage of warmstarts, 
	# this run will be faster
	lambda = 105
	@time solve!(problem, warmstart=true)
	