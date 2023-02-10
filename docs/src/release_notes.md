# Release notes

## v0.15.3 (February 11, 2023)

 * Add support for LDLFactorizations v0.10 [#496](https://github.com/jump-dev/Convex.jl/pull/496).
 * Replace `randn(m, 1)` with `randn(m)` to be more Julian [#498](https://github.com/jump-dev/Convex.jl/pull/498).
 * Add support for indexing expressions with `CartesianIndex` [#500](https://github.com/jump-dev/Convex.jl/pull/500).

## v0.15.2 (August 10, 2022)

 * Add support for LDLFactorizations v0.9 [#493](https://github.com/jump-dev/Convex.jl/pull/493).
 * Fix use of deprecated functions from `AbstractTrees` [#494](https://github.com/jump-dev/Convex.jl/pull/494).

## v0.15.1 (March 28, 2022)
* Use `OrderedDict` internally for reproducible results, issue: [#488](https://github.com/jump-dev/Convex.jl/issues/488), fix: [#489](https://github.com/jump-dev/Convex.jl/pull/489).

## v0.15.0 (March 2, 2022)

### Breaking changes

 * Minimum required version of Julia is now v1.6
 * Updated to MathOptInterface v1.0
   * As a consequence, many previously deprecated solver calls may stop working.
     For example, instead of `() -> SCS.Optimizer(verbose = 0)`, use
     `MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0)`.

## v0.14.18 (November 14, 2021)

* Fix typo in `logisticloss` for length-1 expressions which caused errors (reported in [#458](https://github.com/jump-dev/Convex.jl/issues/458), fixed in [#469](https://github.com/jump-dev/Convex.jl/pull/469)).

## v0.14.17 (November 14, 2021)

* Updated to become compatible with MathOptInterface v0.10, which enables compatibility with the latest version of many solvers
  ([#467](https://github.com/jump-dev/Convex.jl/pull/467), [#468](https://github.com/jump-dev/Convex.jl/pull/468)).

## v0.14.16 (September 25, 2021)

* Improve numerical stability when evaluating `logsumexp` ([#457](https://github.com/jump-dev/Convex.jl/pull/462)). Thanks `@JinraeKim`!

## v0.14.15 (September 15, 2021)

* Use sparse factorization for checking for positive semi-definiteness in `quadform` when possible ([#457](https://github.com/jump-dev/Convex.jl/pull/457)). Thanks `@mtanneau`!
* Add `assume_psd=false` argument to skip checking for positive semi-definiteness in `quadform` ([#456](https://github.com/jump-dev/Convex.jl/pull/456)).


## v0.14.14 (September 8, 2021)

* Increase the tolerance used in checking if a matrix is positive-semi definite in `quadform` ([#453](https://github.com/jump-dev/Convex.jl/pull/453)). Thanks `@numbermaniac`!

## v0.14.13 (July 25, 2021)

* fix `quadform` for positive semi-definite matrices (fixes a regression introduced in v0.14.11 that required strictly positive semi-definite inputs) [#450](https://github.com/jump-dev/Convex.jl/pull/450).

## v0.14.12 (July 19, 2021)

* fix size of result of `evaluate` on `IndexAtom`s [#448](https://github.com/jump-dev/Convex.jl/pull/448). Thanks `@hurak`!

## v0.14.11 (July 5, 2021)

* fix `quadform` in the complex case [#444](https://github.com/jump-dev/Convex.jl/pull/444). Thanks `@lrnv`!

## v0.14.10 (May 20, 2021)

* declare compatibility with BenchmarkTools v1.0 [#441](https://github.com/jump-dev/Convex.jl/pull/441)

## v0.14.9 (May 18, 2021)

* fix some tests in `lp_dual_abs_atom` [#439](https://github.com/jump-dev/Convex.jl/pull/439). Thanks `@moehle`!

## v0.14.8 (May 4, 2021)

* a complete port of [cvxquad](https://github.com/hfawzi/cvxquad) thanks to `@dstahlke`, yielding new functions `quantum_relative_entropy`, `quantum_entropy`, `trace_logm`, `trace_mpower`, and `lieb_ando`, and cones `GeomMeanHypoCone`, `GeomMeanEpiCone`, and `RelativeEntropyEpiCone` ([#418](https://github.com/jump-dev/Convex.jl/pull/418)). Thanks a ton for the awesome contribution `@dstahlke`!

## v0.14.7 (April 22, 2021)

* declare compatibility with BenchmarkTools v0.7 [#434](https://github.com/jump-dev/Convex.jl/pull/434)

## v0.14.6 (March 28, 2021)

* Use `MOI.instantiate` to create the optimizer, which allows users to pass an [`MOI.OptimizerWithAttributes`](https://jump.dev/MathOptInterface.jl/stable/apireference/#MathOptInterface.OptimizerWithAttributes) to configure solver settings [#431](https://github.com/jump-dev/Convex.jl/pull/431). Thanks `@odow`!

## v0.14.5 (March 14, 2021)

* allow `sumlargest(x,k)`, `sumsmallest(x,k)`, and `sumlargesteigs(x,k)` for `k=0` (simply returns `Constant(0)`). ([#429](https://github.com/jump-dev/Convex.jl/pull/429)).

## v0.14.4 (March 14, 2021)

* fixed a bug where the values of variables were being converted to `Float64` even if the problem was solved in high precision. ([#427](https://github.com/jump-dev/Convex.jl/pull/427)).

## v0.14.3 (March 10, 2021)

* update compatibility bounds for BenchmarkTools 0.6

## v0.14.2 (February 15, 2021)

* added lasso, ridge, and elastic net regression examples ([#420](https://github.com/jump-dev/Convex.jl/pull/420)). Thanks to `@PaulSoderlind`!

## v0.14.1 (January 24, 2021)

* there was a bug causing `conj` to act in-place (reported in [#416](https://github.com/jump-dev/Convex.jl/issues/416)), which has been fixed ([#417](https://github.com/jump-dev/Convex.jl/pull/417)). This bug appears to have existed since the introduction of `conj` in Convex.jl v0.5.0.

## v0.14.0 (January 17, 2021)

### Breaking changes

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

### Other changes

* updated `nuclearnorm` and `sumlargesteigs` to allow complex variables, and allow the argument of `sumlargesteigs` to be non-positive-semi-definite ([#409](https://github.com/jump-dev/Convex.jl/pull/409)). Thanks to `@dstahlke`!

## v0.13.8 (December 2, 2020)

* add unary `+` for `Sign` and `ComplexSign` to allow single-argument `hcat` and `vcat` to work ([#405](https://github.com/jump-dev/Convex.jl/pull/405)). Thanks to `@dstahlke`!

## v0.13.7 (September 11, 2020)

* fix [#403](https://github.com/jump-dev/Convex.jl/issues/403) by adding the keyword argument `silent_solver` to `solve!`.

## v0.13.6 (September 8, 2020)

* fix [#401](https://github.com/jump-dev/Convex.jl/issues/401) by allowing `diagm(x)`.

## v0.13.5 (August 25, 2020)

* fix [#398](https://github.com/jump-dev/Convex.jl/issues/398) by allowing `fix!`'d variables in `quadform`.

## v0.13.4 (July 27, 2020)

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

## v0.13.3 (March 22, 2020)

* Make [`add_constraint!`](https://github.com/jump-dev/Convex.jl/pull/381)
  actually add the constraint to the problem.

## v0.13.2 (March 14, 2020)

* Add [`Convex.MAXDIGITS`](https://github.com/jump-dev/Convex.jl/pull/379). Thanks to `@riccardomurri`!

## v0.13.1 (March 6, 2020)

* Allow disabling DCP warnings ([#372](https://github.com/JuliaOpt/Convex.jl/pull/372))
* Restore export of `Constraint` ([#371](https://github.com/JuliaOpt/Convex.jl/pull/371))

## v0.13.0 (February 28, 2020)

### Major changes

* The intermediate layer has changed from MathProgBase.jl to
  [MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)
  ([#330](https://github.com/JuliaOpt/Convex.jl/pull/330)). To solve problems,
  one should pass a MathOptInterface optimizer constructor, such as
  `SCS.Optimizer`, or `MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0)`.
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
