#############################################################################
# sum.jl
# Handles sums of expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.sum
export sum

### Sum Atom
type SumAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function SumAtom(x::AbstractExpr)
    children = (x,)
    return new(:sum, hash(children), children, (1, 1))
  end
end

function sign(x::SumAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::SumAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::SumAtom)
  return ConstVexity()
end

function evaluate(x::SumAtom)
  return sum(evaluate(x.children[1]))
end

# Suppose x was of the form
# x = Ay where A was a coefficient. Then sum(x) can also be considered
# sum(A, 1) * y
function conic_form(x::SumAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    objective, constraints = conic_form(x.children[1], unique_constr)
    new_obj = ConicObj(copy(objective))
    for var in keys(new_obj)
      new_obj[var] = sum(new_obj[var], 1)
    end
    unique_constr[(x.head, x.children_hash)] = new_obj, constraints
  end
  return unique_constr[(x.head, x.children_hash)]
end

sum(x::AbstractExpr) = SumAtom(x)

function sum(x::AbstractExpr, dimension::Int64)
  if dimension == 1
    return Constant(ones(1, x.size[1]), Positive()) * x
  elseif dimension == 2
    return x * Constant(ones(x.size[2], 1), Positive())
  else
    error("Sum not implemented for dimension $dimension")
  end
end
