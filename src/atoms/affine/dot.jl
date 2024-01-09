function ismatrix(x::AbstractExpr)
    return (s = size(x); length(s) == 2 && s[1] > 1 && s[2] > 1)
end
ismatrix(::AbstractMatrix) = true
ismatrix(::Any) = false

# NOTE: Using asvec avoids broadcast-specific behaviors that we want to avoid, such
# as extending singleton dimensions. We need to ensure that the inputs have the same
# length, which broadcast will check for us if both inputs are vectors.
asvec(x) = convert(AbstractExpr, ismatrix(x) ? vec(x) : x)
_vecdot(x, y) = sum(broadcast(*, conj(asvec(x)), asvec(y)))

LinearAlgebra.dot(x::AbstractExpr, y::AbstractExpr) = _vecdot(x, y)
LinearAlgebra.dot(x::Value, y::AbstractExpr) = _vecdot(x, y)
LinearAlgebra.dot(x::AbstractExpr, y::Value) = _vecdot(x, y)
