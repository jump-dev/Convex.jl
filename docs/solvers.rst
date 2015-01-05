=====================================
Solvers
=====================================

By default, Convex.jl uses `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_ to solve linear programs and second-order cone programs (SOCPs).

If you wish to solve semidefinite programs (SDPs) or exponential programs, you will need to install `SCS <https://github.com/JuliaOpt/SCS.jl>`_, which can be done by :code:`Pkg.add("SCS")`.

Any other solver in `JuliaOpt <http://www.juliaopt.org/>`_ may also be used, so long as it supports the conic constraints used to express the problem. Currently, other solvers will be most useful in solving linear programs (LPs) and mixed integer linear programs (MILPs). Mosek and Gurobi can be used to solve SOCPs while Mosek can also solve SDPs. A detailed table can be found `here <http://www.juliaopt.org/>`_ describing which solvers can solve which kind of problems.

Installing these solvers is very simple. Just follow the instructions in each solver's documentation.

To use a specific solver, you can use the following syntax
::

	solve!(p, GurobiSolver())
	solve!(p, MosekSolver())
	solve!(p, GLPKSolverMIP())
	solve!(p, GLPKSolverLP())
	solve!(p, ECOSSolver())
	solve!(p, SCSSolver())

Of course, the solver must be installed first. For example, we can use GLPK to solve a MILP
::

	using GLPKMathProgInterface
	solve!(p, GLPKSolverMIP())

You can set or see the current default solver by
::

	get_default_solver()
	using Gurobi
	set_default_solver(GurobiSolver()) # or set_default_solver(ECOSSolver(verbose=0))
	# Now Gurobi will be used by default as a solver

Many of the solvers such as ECOS or SCS also allow options to be passed in. More details can be found on the specific solver's page. For example, if we wish to increase the maximum number of iterations for ECOS or SCS, we can do so by
::

	using ECOS
	solve!(p, ECOSSolver(maxit=10000))
	using SCS
	solve!(p, SCSSolver(max_iters=10000))
