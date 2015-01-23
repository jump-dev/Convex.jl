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
  else
    return rational_norm(convert(Rational, p))
  end
end

function vecnorm(x::AbstractExpr, p=2)
  vec_x = vec(x)
  return norm(vec_x, p)
end
