# API

!!! info
    See [Supported operations](@ref) for a list of the operations supported by
    Convex.

```@docs
Convex.AbstractVariable
Convex._value
Convex.set_value!
Convex.get_constraints
Convex.add_constraint!
Convex.vexity
Convex.vexity!
Base.sign(x::Convex.AbstractVariable)
Convex.sign!
Convex.VarType
Convex.vartype
Convex.vartype!
Convex.fix!
Convex.free!
Convex.evaluate
Convex.solve!
Convex.MAXDEPTH
Convex.MAXWIDTH
Convex.MAXDIGITS
Convex.ProblemDepot.run_tests
Convex.ProblemDepot.benchmark_suite
Convex.ProblemDepot.foreach_problem
Convex.ProblemDepot.PROBLEMS
Convex.conic_form!
Convex.new_conic_form!
Convex.write_to_file
```
