#############################################################################
# add_subtract.jl
# Handles unary negation, addition and subtraction of variables, constants
# and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export abs
export sign, curvature, monotonicity, evaluate

### Absolute Value

type AbsAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function AbsAtom(x::AbstractExpr)
    children = (x,)
    return new(:abs, hash(children), children, x.size)
  end
end

function sign(x::AbsAtom)
  return Positive()
end

function monotonicity(x::AbsAtom)
  return (Nondecreasing()*sign(x.children[1]),)
end

function curvature(x::AbsAtom)
  return Convex()
end

function evaluate(x::AbsAtom)
  return abs(evaluate(x.children[1]))
end

abs(x::AbstractExpr) = AbsAtom(x)

function dual_conic_form(e::AbsAtom)
  x = e.children[1]
  t = Variable(size(x))
  p = Problem(:minimize, t, x<=t, x>=-t)
  return dual_conic_form(p)
end