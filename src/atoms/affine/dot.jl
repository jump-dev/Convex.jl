import LinearAlgebra.dot
export vecdot, dot


vecdot(x::AbstractExpr, y::AbstractExpr) = sum(broadcast(*, x, y))
vecdot(x::Value, y::AbstractExpr) = sum(broadcast(*, Constant(x), y))
vecdot(x::AbstractExpr, y::Value) = sum(broadcast(*, x, Constant(y)))

dot(x::AbstractExpr, y::AbstractExpr) = (ismatrix(x) || ismatrix(y)) ? error("dot not implemented for matrices. perhaps you're looking for vecdot?") : vecdot(x, y)
dot(x::Value, y::AbstractExpr) = (ismatrix(x) || ismatrix(y)) ? error("dot not implemented for matrices. perhaps you're looking for vecdot?") : vecdot(x, y)
dot(x::AbstractExpr, y::Value) = (ismatrix(x) || ismatrix(y)) ? error("dot not implemented for matrices. perhaps you're looking for vecdot?") : vecdot(x, y)

# tests if an array is a matrix (2D array) with both dimensions of size > 1
function ismatrix(x)
	sz = size(x)
	if length(sz) != 2
		return false
	else
		for s in sz
			if s == 1
				return false
			end
		end
	end
	return true
end
