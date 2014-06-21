import Base.norm, Base.vecnorm
export norm, norm_inf, norm_2, norm_1, norm_fro, vecnorm

function check_size_norm(x::AbstractCvxExpr)
  if x.size[1] > 1 && x.size[2] > 1
    error("norm(x) not supported when x has size $(x.size)")
  end
end

function norm_inf(x::AbstractCvxExpr)
  return max(abs(x))
end

function norm_1(x::AbstractCvxExpr)
  return sum(abs(x))
end

function norm_2(x::AbstractCvxExpr)
  check_size_norm(x)
  vexity = promote_vexity_norm(x)
  this = CvxExpr(:norm_2, [x], vexity, :pos, (1, 1))
  cone_size = get_vectorized_size(x) + 1

  coeffs1 = spzeros(cone_size, 1)
  coeffs1[1] = -1
  coeffs2 = [spzeros(1, get_vectorized_size(x)); -speye(get_vectorized_size(x))]
  coeffs =  VecOrMatOrSparse[coeffs1, coeffs2]
  vars = [this.uid, x.uid]
  constant = zeros(cone_size, 1)

  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, false, true)]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.norm(x.evaluate(), 2)
  return this
end

function norm(x::AbstractCvxExpr, p=2)
  if p == 1
    return norm_1(x)
  elseif p == 2
    return norm_2(x)
  elseif p == Inf
    return norm_inf(x)
  elseif p == :fro
    return norm_2(vec(x))
  else
    error("Norm $p not defined")
  end
end

norm_fro(x::AbstractCvxExpr) = norm_2(vec(x))

function vecnorm(x::AbstractCvxExpr, p=2)
  vec_x = vec(x)
  return norm(vec_x, p)
end
