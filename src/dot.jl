import Base.dot
export dot

dot(x::AbstractExpr, y::AbstractExpr) = x' * y
dot(x::Value, y::AbstractExpr) = Constant(x') * y
dot(x::AbstractExpr, y::Value) = x' * Constant(y)
