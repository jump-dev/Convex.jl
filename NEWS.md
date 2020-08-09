# Changes in v0.14.0

Breaking changes:

* `x + A` will error if `x` is a scalar variable and `A` is an array. Instead, use `x * ones(size(A)) + A`. 
* The objective function can no longer be a constant; e.g., `maximize(1, constraints)` is no longer the same as `satisfy(constraints)`. One needs to use `satisfy` directly for feasability problems.
* The `RelativeEntropyAtom` now returns a scalar value instead of elementwise values. This does not affect the result of `relative_entropy`.
* SDP and exponential constraints now have dual values populated

# Changes in v0.13.4

* You can now create your own variable types by subtyping `AbstractVariable`.
  See the
  [docs](https://www.juliaopt.org/Convex.jl/dev/advanced/#Custom-Variable-Types-1)
  for more information. You can also add constraints directly to a variable
  using `add_constraint!`. ([#358](https://github.com/JuliaOpt/Convex.jl/pull/358))
* Accessors `vexity(x::Variable)`, `sign(x::Variable)`, and
  `evaluate(x::Variable)` should now be the preferred way to access properties
  of a variable; likewise use `set_value!` to set the initial value of a
  variable. ([#358](https://github.com/JuliaOpt/Convex.jl/pull/358))
* To create integer or binary constraints, use the `VarType` enum (e.g.
  `Variable(BinVar)`). Access or set this via `vartype` and `vartype!`.
  ([#358](https://github.com/JuliaOpt/Convex.jl/pull/358))

# Changes in v0.13.3

* Make [`add_constraint!`](https://github.com/jump-dev/Convex.jl/pull/381)
  actually add the constraint to the problem.

# Changes in v0.13.2

* Add [`Convex.MAXDIGITS`](https://github.com/jump-dev/Convex.jl/pull/379)

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
