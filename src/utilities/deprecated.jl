using Base.depwarn

Base.@deprecate norm(x::AbstractExpr, p::Symbol) vecnorm(x, 2)