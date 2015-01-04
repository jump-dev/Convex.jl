=====================================
Quick Tutorial
=====================================

Consider a constrained least squares problem:

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \|Ax - b\|_2^2 \\
    \mbox{subject to} & x >= 0
  \end{array}

where the unknown variable :math:`x` is a vector of size :math:`n`. The values for :math:`A`, :math:`b` are given and have sizes :math:`m\times n`, :math:`m` respectively.

::

	# Let us first make the Convex.jl module available
	Pkg.add('Convex')

	# Generate random problem data
	m = 4;	n = 5
	A = randn(m, n); b = randn(m, 1)

	# Create a (column vector) variable of size n x 1.
	x = Variable(n) # or x = Variable(n, 1)

	# The problem is to minimize ||Ax - b||^2 subject to x >= 0
	# This can be done by: minimize(objective, constraints)
	problem = minimize(sum_squares(A * x + b), [x >= 0])

	# You can add more constraints at any time by
	problem.constraints += [x <= 1, 0.5 <= 2 * x]

	# Solve the problem by calling solve!
	solve!(problem)

	# Check the status of the problem
	problem.status # :Optimal, :Infeasible, :Unbounded etc.

	# Get the optimum value
	problem.optval

.. Get the dual value
.. problem.constraints[1].dual_value

	# Optimal value of variable x or expression sum_squares(A * x + b)
	evaluate(x)
	evaluate(sum_squares(A * x + b))
