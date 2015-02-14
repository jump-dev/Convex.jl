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
