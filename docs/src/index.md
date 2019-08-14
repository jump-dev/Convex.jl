Convex.jl - Convex Optimization in Julia
========================================

Convex.jl is a Julia package for [Disciplined Convex
Programming](http://dcp.stanford.edu/) (DCP). Convex.jl makes it easy to
describe optimization problems in a natural, mathematical syntax, and to
solve those problems using a variety of different (commercial and
open-source) solvers. Convex.jl can solve

> -   linear programs
> -   mixed-integer linear programs and mixed-integer second-order cone
>     programs
> -   dcp-compliant convex programs including
>     -   second-order cone programs (SOCP)
>     -   exponential cone programs
>     -   semidefinite programs (SDP)

Convex.jl supports many solvers, including
[Mosek](https://github.com/JuliaOpt/Mosek.jl),
[Gurobi](https://github.com/JuliaOpt/gurobi.jl),
[ECOS](https://github.com/JuliaOpt/ECOS.jl),
[SCS](https://github.com/karanveerm/SCS.jl) and
[GLPK](https://github.com/JuliaOpt/GLPK.jl), through the
[MathProgBase](http://mathprogbasejl.readthedocs.org/en/latest/)
interface.

Note that Convex.jl was previously called CVX.jl. This package is under
active development; we welcome bug reports and feature requests. For
usage questions, please contact us via the [Julia Discourse](https://discourse.julialang.org/c/domain/opt).
