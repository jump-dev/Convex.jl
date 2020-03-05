# Changes in v0.13.1

* Allow disabling DCP warnings ([#372](https://github.com/JuliaOpt/Convex.jl/pull/372))
* Restore export of `Constraint` ([#371](https://github.com/JuliaOpt/Convex.jl/pull/371))

# Major changes in v0.13.0

* The intermediate layer has changed from MathProgBase.jl to
  [MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)
  ([#330](https://github.com/JuliaOpt/Convex.jl/pull/330)). To solve problems,
  one should pass a MathOptInterface optimizer constructor, such as
  `SCS.Optimizer`, or `() -> SCS.Optimizer(verbose=false)`.
* `lambdamin` and `lambdamax` have been deprecated in favor of `eigmin` and
  `eigmax`. ([#357](https://github.com/JuliaOpt/Convex.jl/pull/357))
* Many "internal" functions and types are no longer exported, such as the atoms,
  types corresponding to constraints and vexities, etc.
  ([#357](https://github.com/JuliaOpt/Convex.jl/pull/357))
* `evaluate(x::Variable)` and `evaluate(c::Constant)` now return scalars and
  vectors as appropriate, instead of `(1,1)`- and `(d,1)`-matrices.
  ([#359](https://github.com/JuliaOpt/Convex.jl/pull/359)). This affects
  functions which used to return `(1,1)`-matrices; e.g., now
  `evaluate(quadform(...))` yields a scalar.
