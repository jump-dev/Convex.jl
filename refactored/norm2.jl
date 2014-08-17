#############################################################################
# norm2.jl
# Handles the euclidean norm (also called frobenius norm for matrices)
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export EucNormAtom, norm2
export sign, monotonicity, curvature, dual_conic_form


type EucNormAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function EucNormAtom(x::AbstractExpr)
    children = (x,)
    return new(:norm2, hash(children), children, (1, 1))
  end

end

function sign(x::EucNormAtom)
  return Positive()
end

function monotonicity(x::EucNormAtom)
  return (sign(x.children[1]) * Nondecreasing(),)
end

function curvature(x::EucNormAtom)
  return Convex()
end

## Create a new variable euc_norm to represent the norm
## Additionally, create the second order conic constraint (euc_norm, x) in SOC
function dual_conic_form(x::EucNormAtom)
  euc_norm = Variable()
  objective, constraints = dual_conic_form(euc_norm)
  child_objective, child_constraints = dual_conic_form(x.children[1])
  append!(constraints, child_constraints)
  soc_constraint = ConicConstr([objective, child_objective], :SOC, [1, get_vectorized_size(x.children[1])])
  push!(constraints, soc_constraint)
  return objective, constraints
end

norm2(x::AbstractExpr) = EucNormAtom(x)

