using Base.depwarn

# Remove underscores

# SOC atoms
Base.@deprecate geo_mean geomean
Base.@deprecate norm_2 norm2
Base.@deprecate quad_form quadform
Base.@deprecate quad_over_lin quadoverlin
Base.@deprecate rational_norm rationalnorm
Base.@deprecate sum_squares sumsquares
Base.@deprecate square_modulus squaremodulus

# SDP atoms
Base.@deprecate lambda_max lambdamax
Base.@deprecate lambda_min lambdamin
Base.@deprecate matrix_frac matrixfrac
Base.@deprecate nuclear_norm nuclearnorm
Base.@deprecate operator_norm operatornorm

# other atoms
Base.@deprecate dot_sort dotsort
Base.@deprecate sum_largest sumlargest
Base.@deprecate sum_smallest sumsmallest

# broadcasting
#if VERSION<v"0.6.0"
#  .*(x::Value, y::AbstractExpr) = broadcast(*, x, y)
#  .*(x::AbstractExpr, y::Value) = broadcast(*, x, y)
#  .*(x::AbstractExpr, y::AbstractExpr) = broadcast(*, x, y)
#  ./(x::Value, y::AbstractExpr) = broadcast(/, x, y)
#  ./(x::AbstractExpr, y::Value) = broadcast(/, x, y)
#  ./(x::AbstractExpr, y::AbstractExpr) = broadcast(/, x, y)
#  .^(x::AbstractExpr, y::Value) = broadcast(^, x, y)
#end
