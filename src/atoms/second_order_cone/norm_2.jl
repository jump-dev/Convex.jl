#############################################################################
# norm_2.jl
# Handles the euclidean norm (also called frobenius norm for matrices)
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.vecnorm
export EucNormAtom, norm_2, vecnorm
export sign, monotonicity, curvature, conic_form


type EucNormAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function EucNormAtom(x::AbstractExpr)
    children = (x,)
    return new(:norm_2, hash(children), children, (1, 1))
  end
end

function sign(x::EucNormAtom)
  return Positive()
end

function monotonicity(x::EucNormAtom)
  return (sign(x.children[1]) * Nondecreasing(),)
end

function curvature(x::EucNormAtom)
  return ConvexVexity()
end

## Create a new variable euc_norm to represent the norm
## Additionally, create the second order conic constraint (euc_norm, x) in SOC
function conic_form(x::EucNormAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    euc_norm = Variable()
    objective = conic_form(euc_norm, unique_conic_forms)
    conic_form(SOCConstraint(euc_norm, x.children[1]), unique_conic_forms)
    add_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

norm_2(x::AbstractExpr) = EucNormAtom(x)
