import Base.norm, Base.vecnorm
export norm_inf, norm, norm_1, vecnorm

norm_inf(x::AbstractExpr) = maximum(abs(x))
norm_1(x::AbstractExpr) = sum(abs(x))
norm_fro(x::AbstractExpr) = norm_2(vec(x))

# behavior of norm should be consistent with julia:
# * vector norms for vectors
# * operator norms for matrices
# * norm(x, :fro) == vecnorm(x, 2)
function norm(x::AbstractExpr, p=2)
  if length(size(x)) <= 1 || minimum(size(x))==1
    # x is a vector
    if p == 1
      return norm_1(x)
    elseif p == 2
      return norm_2(x)
    elseif p == Inf
      return norm_inf(x)
    elseif p > 1
      # TODO: allow tolerance in the rationalize step
      return rational_norm(x, rationalize(Int64, p))
    else
      error("vector p-norms not defined for p < 1")
    end
  else 
    # x is a matrix
    if p == 1
      return sum(x, 1)
    elseif p == 2
      return operator_norm(x)
    elseif p == Inf
      return sum(x, 2)
    elseif p == :fro
      return vecnorm(x, 2)
    else
      error("matrix p-norms only defined for p = 1, 2, and Inf")
    end
  end
end

function vecnorm(x::AbstractExpr, p=2)
  return norm(vec(x), p)
end
