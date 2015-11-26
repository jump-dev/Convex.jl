#############################################################################
# operatornorm.jl
# Handles matrix operator norm (the maximum singular value of a matrix)
# and creates the alias sigmamax
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export operatornorm, sigmamax

### Operator norm

type OperatorNormAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}

  function OperatorNormAtom(x::AbstractExpr)
    children = (x,)
    return new(:operatornorm, hash(children), children, (1,1))
  end
end

function sign(x::OperatorNormAtom)
  return Positive()
end

# The monotonicity
function monotonicity(x::OperatorNormAtom)
  return (NoMonotonicity(),)
end

function curvature(x::OperatorNormAtom)
  return ConvexVexity()
end

# in julia, `norm` on matrices is the operator norm
function evaluate(x::OperatorNormAtom)
  norm(evaluate(x.children[1]), 2)
end

operatornorm(x::AbstractExpr) = OperatorNormAtom(x)
sigmamax(x::AbstractExpr) = OperatorNormAtom(x)

# Create the equivalent conic problem:
#   minimize t
#   subject to
#            [tI_m A; A' tI_n] is positive semidefinite
# see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
# http://arxiv.org/pdf/0706.4138v1.pdf
function conic_form!(x::OperatorNormAtom, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    m, n = size(A)
    t = Variable()
    p = minimize(t, [t*speye(m) A; A' t*speye(n)] ⪰ 0)
    cache_conic_form!(unique_conic_forms, x, p)
  end
  return get_conic_form(unique_conic_forms, x)
end
