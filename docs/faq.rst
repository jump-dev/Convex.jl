=====================================
FAQ
=====================================

.. contents::
  :local:
  :backlinks: none
  :depth: 1

Where can I get help?
--------------------------------------------------
For usage questions, please contact us via the `JuliaOpt mailing list <https://groups.google.com/forum/#!forum/julia-opt>`_.
If you're running into bugs or have feature requests, please use the `Github Issue Tracker <https://github.com/cvxgrp/Convex.jl/issues>`_.

How does Convex.jl differ from JuMP?
------------------------------------
Convex.jl and JuMP are both modelling languages for mathematical programming embedded in Julia, and both
interface with solvers via the MathProgBase interface, so many of the same solvers are available in both.
Convex.jl converts problems to a standard conic form. This approach requires (and certifies) that the problem
is convex and DCP compliant, and guarantees global optimality of the resulting solution.
JuMP allows nonlinear programming through an interface that learns about functions via their derivatives.
This approach is more flexible (for example, you can optimize non-convex functions), but can't
guarantee global optimality if your function is not convex, or warn you if you've entered a non-convex formulation.

For linear programming, the difference is more stylistic. JuMP's syntax is scalar-based and similar to AMPL and GAMS
making it easy and fast to create constraints by indexing and summation (like :code:`sum{x[i], i=1:numLocation}`).
Convex.jl allows (and prioritizes) linear algebraic and functional constructions (like :code:`max(x,y) < A*z`);
indexing and summation are also supported in Convex.jl, but are somewhat slower than in JuMP.
JuMP also lets you efficiently solve a sequence of problems when new constraints are added 
or when coefficients are modified, 
whereas Convex.jl parses the problem again whenever the `solve!` method is called.

Where can I learn more about Convex Optimization?
--------------------------------------------------
See the freely available book `Convex Optimization <http://web.stanford.edu/~boyd/cvxbook/>`_ by Boyd and Vandenberghe for general background on convex optimization.
For help understanding the rules of Disciplined Convex Programming, we recommend this `DCP tutorial website <http://dcp.stanford.edu/>`_.

Are there similar packages available in other languages?
-----------------------------------------------------------
Indeed! You might use `CVXPY <http://www.cvxpy.org>`_ in Python, or `CVX <http://cvxr.com/>`_ in Matlab.

How does Convex.jl work?
-----------------------------------------------------------
For a detailed discussion of how Convex.jl works, see `our paper <http://www.arxiv.org/abs/1410.4821>`_.

How do I cite this package?
---------------------------------------
If you use Convex.jl for published work, we encourage you to cite the software using the following BibTeX citation:
::

	@article{convexjl,
	 title = {Convex Optimization in {J}ulia},
	 author ={Udell, Madeleine and Mohan, Karanveer and Zeng, David and Hong, Jenny and Diamond, Steven and Boyd, Stephen},
	 year = {2014},
	 journal = {SC14 Workshop on High Performance Technical Computing in Dynamic Languages},
	 archivePrefix = "arXiv",
	 eprint = {1410.4821},
	 primaryClass = "math-oc",
	}
