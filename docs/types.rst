=====================================
Basic Types
=====================================

The basic building block of Convex.jl is called an *expression*, which can represent a variable, a constant, or a function of another expression. We discuss each kind of expression in turn.

Variables
=========
The simplest kind of expression in Convex.jl is a variable. Variables in Convex.jl are declared using the `Variable` keyword, along with the dimensions of the variable.
::

	# Scalar variable
	x = Variable()

	# Column vector variable
	x = Variable(5)

	# Matrix variable
	x = Variable(4, 6)

Variables may also be declared as having special properties, such as being

  * (entrywise) positive: :code:`x = Variable(4, Positive())`
  * (entrywise) negative: :code:`x = Variable(4, Negative())`
  * integral: :code:`x = Variable(4, :Int)`
  * binary: :code:`x = Variable(4, :Bin)`
  * (for a matrix) being symmetric, with nonnegative eigenvalues (ie, positive semidefinite): :code:`z = Semidefinite(4)`

Constants
==========
Numbers, vectors, and matrices present in the Julia environment are wrapped automatically into a `Constant` expression when used in a Convex.jl expression.

Expressions
============
Expressions in Convex.jl are formed by applying any *atom* (mathematical function defined in Convex.jl) to variables, constants, and other expressions. For a list of these functions, see `operations <operations.html>`_.
Atoms are applied to expressions using operator overloading. For example, :code:`2+2` calls Julia's built-in addition operator, while :code:`2+x` calls the Convex.jl addition method and returns a Convex.jl expression. Many of the useful language features in Julia, such as arithmetic, array indexing, and matrix transpose are overloaded in Convex.jl so they may be used with variables and expressions just as they are used with native Julia types.

Expressions that are created must be DCP-compliant.
More information on DCP can be found `here <http://dcp.stanford.edu/>`_.
::

	x = Variable(5)
	# The following are all expressions
	y = sum(x)
	z = 4 * x + y
	z_1 = z[1]

Convex.jl allows the values of the expressions to be evaluated directly.
::

	x = Variable()
	y = Variable()
	z = Variable()
	expr = x + y + z
	problem = minimize(expr, x >= 1, y >= x, 4 * z >= y)
	solve!(problem, SCSSolver())

	# Once the problem is solved, we can call evaluate() on expr:
	evaluate(expr)


Constraints
============
*Constraints* in Convex.jl are declared using the standard comparison operators :code:`<=`, :code:`>=`, and :code:`==`.  They specify relations that must hold between two expressions.  Convex.jl does not distinguish between strict and non-strict inequality constraints.
::

	x = Variable(5, 5)
	# Equality constraint
	constraint = x == 0
	# Inequality constraint
	constraint = x >= 1

Note that constraints apply elementwise automatically; that is, :code:`x >= 1` means that :code:`x[i] >= 1` for :code:`i=1:5`. Similarly, if :code:`y = rand(5)`, then :code:`x >= y` means that :code:`x[i] >= 1` for :code:`i=1:5`. In particular, broadcast should not generally be used to constain vectors, i.e., use :code:`x >= y` instead of :code:`x .>= y`. The latter yields a vector of constraints instead of one constraint on the vector.

Matrices can also be constrained to be positive semidefinite.
::

	x = Variable(3, 3)
	y = Variable(3, 1)
	z = Variable()
	# constrain [x y; y' z] to be positive semidefinite
	constraint = ([x y; y' z] in :SDP)
	# or equivalently,
	constraint = ([x y; y' z] ⪰ 0)

Objective
=========
The objective of the problem is a scalar expression to be maximized or minimized by using :code:`maximize` or :code:`minimize` respectively. Feasibility problems can be expressed by either giving a constant as the objective, or using :code:`problem = satisfy(constraints)`.

Problem
========
A *problem* in Convex.jl consists of a *sense* (minimize, maximize, or satisfy), an *objective* (an expression to which the sense verb is to be
applied), and zero or more *constraints* that must be satisfied at the solution.
Problems may be constructed as
::

	problem = minimize(objective, constraints)
	# or
	problem = maximize(objective, constraints)
	# or
	problem = satisfy(constraints)

Constraints can be added at any time before the problem is solved.
::

	# No constraints given
	problem = minimize(objective)
	# Add some constraint
	problem.constraints += constraint
	# Add many more constraints
	problem.constraints += [constraint1, constraint2, ...]

A problem can be solved by calling :code:`solve!`:
::

	solve!(problem, solver)

passing a solver such as :code:`SCSSolver()` from the package :code:`SCS` as the second argument.
After the problem is solved, :code:`problem.status` records the status returned by the optimization solver, and can be :code:`:Optimal`, :code:`:Infeasible`, :code:`:Unbounded`, :code:`:Indeterminate` or :code:`:Error`.
If the status is :code:`:Optimal`, :code:`problem.optval` will record the optimum value of the problem.
The optimal value for each variable :code:`x` participating in the problem can be found in :code:`x.value`.
The optimal value of an expression can be found by calling the :code:`evaluate()` function on the expression as follows: :code:`evaluate(expr)`.

.. The dual values are stored with the respective constraints and can be accessed as :code:`problem.constraints[idx].dual_value`.
