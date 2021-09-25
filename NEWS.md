# Changes in v0.14.16

* Improve numerical stability when evaluating `logsumexp` ([#457](https://github.com/jump-dev/Convex.jl/pull/462)). Thanks @JinraeKim!

# Changes in v0.14.15

* Use sparse factorization for checking for positive semi-definiteness in `quadform` when possible ([#457](https://github.com/jump-dev/Convex.jl/pull/457)). Thanks @mtanneau!
* Add `assume_psd=false` argument to skip checking for positive semi-definiteness in `quadform` ([#456](https://github.com/jump-dev/Convex.jl/pull/456)).


# Changes in v0.14.14

* Increase the tolerance used in checking if a matrix is positive-semi definite in `quadform` ([#453](https://github.com/jump-dev/Convex.jl/pull/453)). Thanks @numbermaniac!

# Changes in v0.14.13

* fix `quadform` for positive semi-definite matrices (fixes a regression introduced in v0.14.11 that required strictly positive semi-definite inputs) [#450](https://github.com/jump-dev/Convex.jl/pull/450).

# Changes in v0.14.12

* fix size of result of `evaluate` on `IndexAtom`s [#448](https://github.com/jump-dev/Convex.jl/pull/448). Thanks @hurak!

# Changes in v0.14.11

* fix `quadform` in the complex case [#444](https://github.com/jump-dev/Convex.jl/pull/444). Thanks @lrnv!

# Changes in v0.14.10

* declare compatibility with BenchmarkTools v1.0 [#441](https://github.com/jump-dev/Convex.jl/pull/441)

# Changes in v0.14.9

* fix some tests in `lp_dual_abs_atom` [#439](https://github.com/jump-dev/Convex.jl/pull/439). Thanks @moehle!

# Changes in v0.14.8

* a complete port of [cvxquad](https://github.com/hfawzi/cvxquad) thanks to @dstahlke, yielding new functions `quantum_relative_entropy`, `quantum_entropy`, `trace_logm`, `trace_mpower`, and `lieb_ando`, and cones `GeomMeanHypoCone`, `GeomMeanEpiCone`, and `RelativeEntropyEpiCone` ([#418](https://github.com/jump-dev/Convex.jl/pull/418)). Thanks a ton for the awesome contribution @dstahlke!

# Changes in v0.14.7

* declare compatibility with BenchmarkTools v0.7 [#434](https://github.com/jump-dev/Convex.jl/pull/434)

# Changes in v0.14.6

* Use `MOI.instantiate` to create the optimizer, which allows users to pass an [`MOI.OptimizerWithAttributes`](https://jump.dev/MathOptInterface.jl/stable/apireference/#MathOptInterface.OptimizerWithAttributes) to configure solver settings [#431](https://github.com/jump-dev/Convex.jl/pull/431). Thanks @odow!

# Changes in v0.14.5

* allow `sumlargest(x,k)`, `sumsmallest(x,k)`, and `sumlargesteigs(x,k)` for `k=0` (simply returns `Constant(0)`). ([#429](https://github.com/jump-dev/Convex.jl/pull/429)).

# Changes in v0.14.4

* fixed a bug where the values of variables were being converted to `Float64` even if the problem was solved in high precision. ([#427](https://github.com/jump-dev/Convex.jl/pull/427)).

# Changes in v0.14.3

* update compatibility bounds for BenchmarkTools 0.6

# Changes in v0.14.2

* added lasso, ridge, and elastic net regression examples ([#420](https://github.com/jump-dev/Convex.jl/pull/420)). Thanks to @PaulSoderlind!

# Changes in v0.14.1

* there was a bug causing `conj` to act in-place (reported in [#416](https://github.com/jump-dev/Convex.jl/issues/416)), which has been fixed ([#417](https://github.com/jump-dev/Convex.jl/pull/417)). This bug appears to have existed since the introduction of `conj` in Convex.jl v0.5.0.

# Changes in v0.14.0

## Breaking changes

* Changes to the `sign` of atoms:
    * The sign of `sumlargesteigs` has been changed from `Positive()` to  `NoSign()`, to allow non-positive-semidefinite inputs ([#409](https://github.com/jump-dev/Convex.jl/pull/409)). This has the potential to break code that required that sign to be positive. If you run into this problem, please file an issue so we can figure out a workaround.
    * The sign of `eigmin` and `eigmax` has been changed from `Positive()` to  `NoSign()` ([#413](https://github.com/jump-dev/Convex.jl/pull/413)). This is a bugfix because in general `eigmin` and `eigmax` do not need to return a positive quantity (for non-positive-semidefinite inputs). Again, this has the potential to break code that required that sign to be positive. If you run into this problem, please file an issue so we can figure out a workaround.
* Removal of deprecations:
    * `lambdamin` and `lambdamax` has been deprecated to `eigmin` and `eigmax` since Convex v0.13.0. This deprecation has been removed, so your code must be updated to call `eigmin` or `eigmax` instead ([#412](https://github.com/jump-dev/Convex.jl/pull/412)).
    * `norm(x, p)` where `x` is a matrix expression has been deprecated to `opnorm(x,p)` since Convex v0.8.0. This deprecation has been removed, so your code must be updated to call `opnorm(x, p)` instead ([#412](https://github.com/jump-dev/Convex.jl/pull/412)). Currently, `norm(x,p)` for a matrix
    expression `x` will error, but in Convex.jl v0.15.0 it will return `norm(vec(x), p)`.
    * `Convex.clearmemory()` has been deprecated and unnecessary since Convex v0.12.5. This deprecation has been removed, so if this function is in your code, just delete it ([#412](https://github.com/jump-dev/Convex.jl/pull/412)).
    * `vecnorm(x, p)` has been deprecated to `norm(vec(x), p)` since Convex v0.8.0. This deprecation has been removed, so your code must be updated to call `norm(vec(x),p)` instead ([#412](https://github.com/jump-dev/Convex.jl/pull/412)).
* Other changes:
    * `Convex.DCP_WARNINGS` was introduced in Convex v0.13.1 to allow turning off Convex.jl's DCP warnings. This has been removed in favor of the function `Convex.emit_dcp_warnings()` ([Commit 481fa02](https://github.com/jump-dev/Convex.jl/commit/481fa02b84bfec6bf7c809ea93d6ba8004193b83)).

## Other changes

* updated `nuclearnorm` and `sumlargesteigs` to allow complex variables, and allow the argument of `sumlargesteigs` to be non-positive-semi-definite ([#409](https://github.com/jump-dev/Convex.jl/pull/409)). Thanks to @dstahlke!

# Changes in v0.13.8

* add unary `+` for `Sign` and `ComplexSign` to allow single-argument `hcat` and `vcat` to work ([#405](https://github.com/jump-dev/Convex.jl/pull/405)). Thanks to @dstahlke!

# Changes in v0.13.7

* fix [#403](https://github.com/jump-dev/Convex.jl/issues/403) by adding the keyword argument `silent_solver` to `solve!`.

# Changes in v0.13.6

* fix [#401](https://github.com/jump-dev/Convex.jl/issues/401) by allowing `diagm(x)`.

# Changes in v0.13.5

* fix [#398](https://github.com/jump-dev/Convex.jl/issues/398) by allowing `fix!`'d variables in `quadform`.

# Changes in v0.13.4

* You can now create your own variable types by subtyping `AbstractVariable`.
  See the
  [docs](https://www.juliaopt.org/Convex.jl/dev/advanced/#Custom-Variable-Types-1)
  for more information. You can also add constraints directly to a variable
  using `add_constraint!` ([#358](https://github.com/JuliaOpt/Convex.jl/pull/358)).
* Accessors `vexity(x::Variable)`, `sign(x::Variable)`, and
  `evaluate(x::Variable)` should now be the preferred way to access properties
  of a variable; likewise use `set_value!` to set the initial value of a
  variable ([#358](https://github.com/JuliaOpt/Convex.jl/pull/358)).
* To create integer or binary constraints, use the `VarType` enum (e.g.
  `Variable(BinVar)`). Access or set this via `vartype` and `vartype!` ([#358](https://github.com/JuliaOpt/Convex.jl/pull/358)).

# Changes in v0.13.3

* Make [`add_constraint!`](https://github.com/jump-dev/Convex.jl/pull/381)
  actually add the constraint to the problem.

# Changes in v0.13.2

* Add [`Convex.MAXDIGITS`](https://github.com/jump-dev/Convex.jl/pull/379). Thanks to @riccardomurri!

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
  `eigmax` ([#357](https://github.com/JuliaOpt/Convex.jl/pull/357)).
* Many "internal" functions and types are no longer exported, such as the atoms,
  types corresponding to constraints and vexities, etc.
  ([#357](https://github.com/JuliaOpt/Convex.jl/pull/357)).
* `evaluate(x::Variable)` and `evaluate(c::Constant)` now return scalars and
  vectors as appropriate, instead of `(1,1)`- and `(d,1)`-matrices
  ([#359](https://github.com/JuliaOpt/Convex.jl/pull/359)). This affects
  functions which used to return `(1,1)`-matrices; e.g., now
  `evaluate(quadform(...))` yields a scalar.
