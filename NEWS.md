# Changes in v0.13.0

* The intermediate layer has changed from MathProgBase.jl to
  [MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)
  ([#330](https://github.com/JuliaOpt/Convex.jl/pull/330)).
* `lambdamin` and `lambdamax` have been deprecated in favor of `eigmin` and
  `eigmax`. ([#357](https://github.com/JuliaOpt/Convex.jl/pull/357))
* `evaluate(x::Variable)` and `evaluate(c::Constant)` now return scalars and
  vectors as appropriate, instead of `(1,1)`- and `(d,1)`-matrices.
  ([#359](https://github.com/JuliaOpt/Convex.jl/pull/359)). This affects
  functions which used to return `(1,1)`-matrices; e.g., now
  `evaluate(quadform(...))` yields a scalar.
