=====================================
Installation
=====================================


Installing Convex.jl is a one step process. Open up Julia and type
::

	Pkg.update()
	Pkg.add("Convex")

This also installs `SCS <https://github.com/JuliaOpt/SCS.jl>`_ as a solver. To solve certain problems such as mixed integer programming problems you will need to install another solver as well, such as `GLPK <https://github.com/JuliaOpt/GLPKMathProgInterface.jl>`_. If you wish to use other solvers, please read the section on `solvers <solvers.html>`_.
