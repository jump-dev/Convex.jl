=====================================
Installation
=====================================


Installing Convex.jl is a one step process. Open up Julia and type
::

	Pkg.update()
	Pkg.add("Convex")

This also installs `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_ as a solver. If you wish to use other solvers, please read the section on `solvers <solvers.html>`_.
