=====================================
Installation
=====================================


Installing Convex.jl is a one step process. Open up Julia and type
::

	Pkg.update()
	Pkg.add("Convex")

This also installs `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_ as a solver. To use certain functions (exponential functions like :code:`exp` and :code:`log`, and semidefinite functions like :code:`operator_norm` and :code:`logdet`) you will need to install another solver as well, such as `SCS <https://github.com/JuliaOpt/SCS.jl>`_. If you wish to use other solvers, please read the section on `solvers <solvers.html>`_.
