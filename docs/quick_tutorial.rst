=====================================
Quick Tutorial
=====================================

Consider a constrained least squares problem

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \|Ax - b\|_2^2 \\
    \mbox{subject to} & x \geq 0
  \end{array}

with variable :math:`x\in \mathbf{R}^{n}`, 
and problem data :math:`A \in \mathbf{R}^{m \times n}`, :math:`b \in \mathbf{R}^{m}`.

This problem can be solved in Convex.jl as follows:
::

	# Make the Convex.jl module available
	using Convex

	# Generate random problem data
	m = 4;	n = 5
	A = randn(m, n); b = randn(m, 1)

	# Create a (column vector) variable of size n x 1.
	x = Variable(n)

	# The problem is to minimize ||Ax - b||^2 subject to x >= 0
	# This can be done by: minimize(objective, constraints)
	problem = minimize(sum_squares(A * x - b), [x >= 0])

	# Solve the problem by calling solve!
	solve!(problem)

	# Check the status of the problem
	problem.status # :Optimal, :Infeasible, :Unbounded etc.

	# Get the optimum value
	problem.optval

.. Get the dual value
.. problem.constraints[1].dual_value

	# Optimal value of variable x or expression sum_squares(A * x - b)
	evaluate(x)
	evaluate(sum_squares(A * x - b))
