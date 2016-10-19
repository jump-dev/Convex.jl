import Base.norm, Base.vecnorm
export norm_inf, norm, norm_1, vecnorm

# deprecate these soon
norm_inf(x::AbstractExpr) = maximum(abs(x))
norm_1(x::AbstractExpr) = sum(abs(x))
norm_fro(x::AbstractExpr) = norm2(vec(x))

# behavior of norm should be consistent with julia:
# * vector norms for vectors
# * operator norms for matrices
function norm(x::AbstractExpr, p::Number=2)
  if length(size(x)) <= 1 || minimum(size(x))==1
    # x is a vector
    if p == 1
      return norm_1(x)
    elseif p == 2
      return norm2(x)
    elseif p == Inf
      return norm_inf(x)
    elseif p > 1
      # TODO: allow tolerance in the rationalize step
      return rationalnorm(x, rationalize(Int, float(p)))
    else
      error("vector p-norms not defined for p < 1")
    end
  else
    # x is a matrix
    if p == 1
      return maximum(sum(abs(x), 1))
    elseif p == 2
      return operatornorm(x)
    elseif p == Inf
      return maximum(sum(abs(x), 2))
    else
      error("matrix p-norms only defined for p = 1, 2, and Inf")
    end
  end
end

function vecnorm(x::AbstractExpr, p::Real=2)
  return norm(vec(x), p)
end
