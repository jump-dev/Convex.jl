# Solvers

Convex.jl transforms each problem into an equivalent [cone
program](http://mathprogbasejl.readthedocs.org/en/latest/conic.html) in
order to pass the problem to a specialized solver. Depending on the
types of functions used in the problem, the conic constraints may
include linear, second-order, exponential, or semidefinite constraints,
as well as any binary or integer constraints placed on the variables.

By default, Convex.jl does not install any solvers. Many users use the solver
[SCS](https://github.com/JuliaOpt/SCS.jl), which is able to solve problems with
linear, second-order cone constraints (SOCPs), exponential constraints and
semidefinite constraints (SDPs). Likewise,
[COSMO](https://github.com/oxfordcontrol/COSMO.jl) is a pure-Julia solver which
can handle every cone that Convex.jl itself supports. Any other solver in
[JuliaOpt](http://www.juliaopt.org/) may also be used, so long as it supports
the conic constraints used to represent the problem. Many other solvers in the
JuliaOpt ecosystem can be used to solve (mixed integer) linear programs (LPs and
MILPs). Mosek and Gurobi can be used to solve SOCPs (even with binary or integer
constraints), and Mosek can also solve SDPs. For up-to-date information about
solver capabilities, please see the table [here](http://www.juliaopt.org/)
describing which solvers can solve which kind of problems. See also
[ConvexTests.jl](https://ericphanson.github.io/ConvexTests.jl/dev/) to see the
results of running test problems with Convex.jl for many solvers.

Installing these solvers is very simple. Just follow the instructions in
the documentation for that solver.

To use a specific solver, you can use the following syntax

```julia
solve!(p, Gurobi.Optimizer)
solve!(p, Mosek.Optimizer)
solve!(p, GLPK.Optimizer)
solve!(p, ECOS.Optimizer)
solve!(p, SCS.Optimizer)
```

(Of course, the solver must be installed first.) For example, we can use
GLPK to solve a MILP:
```julia
using GLPK
solve!(p, GLPK.Optimizer)
```

Many of the solvers also allow options to be passed using
`MOI.OptimizerWithAttributes`. For example, to set a time limit (in
milliseconds) with GLPK, use:
```julia
using Convex, GLPK
const MOI = Convex.MOI

solve!(
    p,
    MOI.OptimizerWithAttributes(GLPK.Optimizer, "tm_lim" => 60_000.0)
)
```

As another example, if we wish to turn off printing for the SCS solver
(that is, run in quiet mode), we can do so as follows:
```julia
using Convex, SCS
const MOI = Convex.MOI

opt = MOI.OptimizerWithAttributes(SCS.Optimizer, MOI.Silent() => false)
solve!(p, opt)
```

Another option is to use the solver-independent `silent_solver`
keyword argument to `solve!`:
```julia
solve!(p, SCS.Optimizer; silent_solver=true)
```

See each solver's documentation for more information on solver-dependent
options.

To turn off the problem status warning issued by Convex when a solver is
not able to solve a problem to optimality, use the keyword argument
`verbose=false` of the solve method:

```julia
solve!(p, SCS.Optimizer, verbose=false)
```
