# Solvers

Convex.jl depends on third-party solvers to solve optimization problems.
Therefore, you will need to install one before you can solve problems with
Convex.jl.

Install a solver using the Julia package manager, replacing `"SCS"` by Julia
package name as appropriate:
```julia
import Pkg
Pkg.add("SCS")
```

The [JuMP documentation](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
has a list of support solvers and a list of problem classes they support.

To use a specific solver, you can use the following syntax:
```@repl solvers
using Convex, SCS
x = Variable();
p = minimize(x, [x >= 1]);
solve!(p, SCS.Optimizer)
```

A different solver can be used by replacing `SCS` as appropriate. For example,
GLPK is a mixed-inter linear solver:
```@repl solvers
using GLPK
solve!(p, GLPK.Optimizer)
```

Many of the solvers also allow options to be passed using
`MOI.OptimizerWithAttributes`. For example, to set a time limit (in
milliseconds) with GLPK, use:
```@repl solvers
import Convex: MOI
solver = MOI.OptimizerWithAttributes(GLPK.Optimizer, "tm_lim" => 60_000.0)
solve!(p, solver)
```

As another example, if we wish to turn off printing for the SCS solver
(that is, run in quiet mode), we can do so as follows:
```@repl solvers
silent_scs = MOI.OptimizerWithAttributes(SCS.Optimizer, MOI.Silent() => true)
solve!(p, silent_scs)
```

Another option is to use the solver-independent `silent_solver` keyword
argument to `solve!`:
```@repl solvers
solve!(p, SCS.Optimizer; silent_solver=true)
```

See each solver's documentation for more information on solver-dependent
options.
