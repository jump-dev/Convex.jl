Installation
============

Installing Convex.jl is a one step process. Open up Julia and type :

```julia
using Pkg
Pkg.update()
Pkg.add("Convex")
```

This does not install any solvers. If you don't have a solver installed
already, you will want to install a solver such as
[SCS](https://github.com/JuliaOpt/SCS.jl) by running :

```julia
Pkg.add("SCS")
```

To solve certain problems such as mixed integer programming problems you
will need to install another solver as well, such as
[GLPK](https://github.com/JuliaOpt/GLPK.jl). If you
wish to use other solvers, please read the section on
[Solvers](@ref).
