=====================================
Solvers
=====================================

Convex.jl transforms each problem into an equivalent `cone program <http://mathprogbasejl.readthedocs.org/en/latest/conic.html>`_ in order to pass the problem to a specialized solver.
Depending on the types of functions used in the problem, the conic constraints may include linear, second-order, exponential, or semidefinite constraints, as well as any binary or integer constraints placed on the variables.

By default, Convex.jl uses the solver `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_, which is able to solve problems with linear and second-order cone constraints (SOCPs).
If you wish to solve problems with semidefinite or exponential cone constraints, you will need to install `SCS <https://github.com/JuliaOpt/SCS.jl>`_.
Any other solver in `JuliaOpt <http://www.juliaopt.org/>`_ may also be used, so long as it supports the conic constraints used to represent the problem.
Most solvers other solvers in the JuliaOpt ecosystem can be used to solve (mixed integer) linear programs (LPs and MILPs).
Mosek and Gurobi can be used to solve SOCPs (even with binary or integer constraints), and Mosek can also solve SDPs.
For up-to-date information about solver capabilities, please see the table `here <http://www.juliaopt.org/>`_
describing which solvers can solve which kind of problems.

Installing these solvers is very simple. Just follow the instructions in the documentation for that solver.

To use a specific solver, you can use the following syntax
::

	solve!(p, GurobiSolver())
	solve!(p, MosekSolver())
	solve!(p, GLPKSolverMIP())
	solve!(p, GLPKSolverLP())
	solve!(p, ECOSSolver())
	solve!(p, SCSSolver())

(Of course, the solver must be installed first.) For example, we can use GLPK to solve a MILP
::

	using GLPKMathProgInterface
	solve!(p, GLPKSolverMIP())

You can set or see the current default solver by
::

	get_default_solver()
	using Gurobi
	set_default_solver(GurobiSolver()) # or set_default_solver(ECOSSolver(verbose=0))
	# Now Gurobi will be used by default as a solver

Many of the solvers also allow options to be passed in. More details can be found in each solver's documentation.
For example, if we wish to increase the maximum number of iterations for ECOS or SCS, we can do so by
::

	using ECOS
	solve!(p, ECOSSolver(maxit=10000))
	using SCS
	solve!(p, SCSSolver(max_iters=10000))
