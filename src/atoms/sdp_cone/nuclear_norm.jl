#############################################################################
# nuclear_norm.jl
# Handles nuclear norm (the sum of the singular values of a matrix),
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export nuclear_norm

### Nuclear norm

type NuclearNormAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::@compat Tuple{AbstractExpr}
  size::@compat Tuple{Int, Int}

  function NuclearNormAtom(x::AbstractExpr)
    children = (x,)
    return new(:nuclear_norm, hash(children), children, (1,1))
  end
end

function sign(x::NuclearNormAtom)
  return Positive()
end

# The monotonicity
function monotonicity(x::NuclearNormAtom)
  return (NoMonotonicity(),)
end

function curvature(x::NuclearNormAtom)
  return ConvexVexity()
end

function evaluate(x::NuclearNormAtom)
  return sum(svdvals(evaluate(x.children[1])))
end

nuclear_norm(x::AbstractExpr) = NuclearNormAtom(x)

# Create the equivalent conic problem:
#   minimize (trace(U) + trace(V))/2
#   subject to
#            [U A; A' V] is positive semidefinite
# see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
# http://arxiv.org/pdf/0706.4138v1.pdf
function conic_form!(x::NuclearNormAtom, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    m, n = size(A)
    U = Variable(m,m)
    V = Variable(n,n)
    p = minimize(.5*(trace(U) + trace(V)), [U A; A' V] ⪰ 0)
    cache_conic_form!(unique_conic_forms, x, p)
  end
  return get_conic_form(unique_conic_forms, x)
end
