import Base.dot
export dot

dot(x::AbstractExpr, y::AbstractExpr) = sum(x .* y)
dot(x::Value, y::AbstractExpr) = sum(Constant(x) .* y)
dot(x::AbstractExpr, y::Value) = sum(x .* Constant(y))
