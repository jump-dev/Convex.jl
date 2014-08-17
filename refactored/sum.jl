#############################################################################
# sum.jl
# Handles sums of expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export sum

### Sums

type SumAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function SumAtom(x::AbstractExpr)
    children = (x,)
    return new(:sum, hash(children), children, x.size)
  end
end

function sign(x::SumAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::SumAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x); this represents the first term
function curvature(x::SumAtom)
  return ConstVexity()
end

function evaluate(x::NegateAtom)
  return sum(evaluate(x.children[1]))
end

sum(x::AbstractExpr) = SumAtom(x)

function dual_conic_form(e::SumAtom)
  objective, constraints = dual_conic_form(e.children[1])
  new_obj = ConicObj(copy(objective.vars_to_coeffs))
  for var in keys(new_obj.vars_to_coeffs)
    new_obj.vars_to_coeffs[var] = sum(new_obj.vars_to_coeffs[var], 1)
  end
  return new_obj, constraints
end