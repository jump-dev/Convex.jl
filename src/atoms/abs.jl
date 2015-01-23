#############################################################################
# abs.jl
# Absolute value of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.abs
export abs
export sign, curvature, monotonicity, evaluate

### Absolute Value

type AbsAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int, Int)

  function AbsAtom(x::AbstractExpr)
    children = (x,)
    return new(:abs, hash(children), children, x.size)
  end
end

function sign(x::AbsAtom)
  return Positive()
end

function monotonicity(x::AbsAtom)
  return (Nondecreasing() * sign(x.children[1]),)
end

function curvature(x::AbsAtom)
  return ConvexVexity()
end

function evaluate(x::AbsAtom)
  return abs(evaluate(x.children[1]))
end

abs(x::AbstractExpr) = AbsAtom(x)

function conic_form!(x::AbsAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    c = x.children[1]
    t = Variable(size(c))
    objective = conic_form!(t, unique_conic_forms)
    conic_form!(c<=t, unique_conic_forms)
    conic_form!(c>=-t, unique_conic_forms)
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end
