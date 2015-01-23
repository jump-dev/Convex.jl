=====================================
Installation
=====================================


Installing Convex.jl is a one step process. Open up Julia and type
::

	Pkg.update()
	Pkg.add("Convex")

This does not install any solvers. If you don't have a solver installed already, you will want to install a solver such as `SCS <https://github.com/JuliaOpt/SCS.jl>`_ by running
::

	Pkg.add("SCS")

To solve certain problems such as mixed integer programming problems you will need to install another solver as well, such as `GLPK <https://github.com/JuliaOpt/GLPKMathProgInterface.jl>`_. If you wish to use other solvers, please read the section on `solvers <solvers.html>`_.
