# Using _vec avoids broadcast-specific behaviors that we want to avoid,
# such as extending singleton dimensions. We need to ensure that the inputs have
# the same length, which broadcast will check for us if both inputs are vectors.
_vec(x::AbstractExpr) = size(x, 1) > 1 && size(x, 2) > 1 ? vec(x) : x

_vec(x::AbstractMatrix) = vec(x)

_vec(x::Any) = x

_expr_vec(x) = convert(AbstractExpr, _vec(x))

_dot(x, y) = sum(broadcast(*, conj(_expr_vec(x)), _expr_vec(y)))

LinearAlgebra.dot(x::AbstractExpr, y::AbstractExpr) = _dot(x, y)

LinearAlgebra.dot(x::Value, y::AbstractExpr) = _dot(x, y)

LinearAlgebra.dot(x::AbstractExpr, y::Value) = _dot(x, y)
