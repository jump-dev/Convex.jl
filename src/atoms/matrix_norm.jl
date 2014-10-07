#############################################################################
# matrix_norm.jl
# Handles nuclear norm (the sum of the singular values of a matrix),
# matrix operator norm (the maximum singular value of a matrix),
# and the trace of a matrix
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.trace
export nuclear_norm, operator_norm, sigma_max, trace

### Nuclear norm

type NuclearNormAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function NuclearNormAtom(x::AbstractExpr)
    children = (x,)
    return new(:nuclear_norm, hash(children), children, x.size)
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
function conic_form(x::NuclearNormAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    A = x.children[1]
    m, n = size(A)
    U = Variable(m,m)
    V = Variable(n,n)
    constraint = isposdef([U A; A' V])
    objective = .5*(trace(U) + trace(V))
    p = minimize(objective, constraint)
    unique_constr[(x.head, x.children_hash)] = conic_form(p)
  end
  return safe_copy(unique_constr[(x.head, x.children_hash)])
end

### Trace

function trace(e::AbstractExpr)
  return sum(diag(e))
end

### Operator norm

type OperatorNormAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function OperatorNormAtom(x::AbstractExpr)
    children = (x,)
    return new(:operator_norm, hash(children), children, x.size)
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

# XXX verify this returns all the eigenvalues even in new versions of julia (>=3.0)
function evaluate(x::OperatorNormAtom)
  svdvals(evaluate(x.children[1]))[1]
end

operator_norm(x::AbstractExpr) = OperatorNormAtom(x)
sigma_max(x::AbstractExpr) = OperatorNormAtom(x)

# Create the equivalent conic problem:
#   minimize t
#   subject to
#            [tI_m A; A' tI_n] is positive semidefinite
# see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
# http://arxiv.org/pdf/0706.4138v1.pdf
function conic_form(x::OperatorNormAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    A = x.children[1]
    m, n = size(A)
    t = Variable()
    constraint = isposdef([t*speye(m) A; A' t*speye(n)])
    objective = t
    p = minimize(objective, constraint)
    unique_constr[(x.head, x.children_hash)] = conic_form(p)
  end
  return safe_copy(unique_constr[(x.head, x.children_hash)])
end
