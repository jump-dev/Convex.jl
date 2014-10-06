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
  return ConvexVexity()
end

function evaluate(x::AbsAtom)
  return abs(evaluate(x.children[1]))
end

abs(x::AbstractExpr) = AbsAtom(x)

function conic_form(x::AbsAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    c = x.children[1]
    t = Variable(size(c))
    objective, constraints = conic_form(t, unique_constr)
    append!(constraints, conic_form(c<=t, unique_constr)[2])
    append!(constraints, conic_form(c>=-t, unique_constr)[2])
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end
