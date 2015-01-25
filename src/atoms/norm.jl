import Base.norm, Base.vecnorm
export norm_inf, norm, norm_1, vecnorm

norm_inf(x::AbstractExpr) = maximum(abs(x))
norm_1(x::AbstractExpr) = sum(abs(x))
norm_fro(x::AbstractExpr) = norm_2(vec(x))

function norm(x::AbstractExpr, p=2)
  if p == 1
    return norm_1(x)
  elseif p == 2
    return norm_2(x)
  elseif p == Inf
    return norm_inf(x)
  elseif p == :fro
    return norm_2(vec(x))
  elseif p > 1
    return rational_norm(x, convert(Rational, p))
  else
    error("p-norms not defined for p < 1")
  end
end

function vecnorm(x::AbstractExpr, p=2)
  vec_x = vec(x)
  return norm(vec_x, p)
end
