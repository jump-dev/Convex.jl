#############################################################################
# minimum.jl
# Compute the minimum value of an array.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.minimum
export minimum

### Minimum Atom
type MinimumAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function MinimumAtom(x::AbstractExpr)
    children = (x,)
    return new(:minimum, hash(children), children, (1, 1))
  end
end

function sign(x::MinimumAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::MinimumAtom)
  return (Nonincreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MinimumAtom)
  return ConcaveVexity()
end

function evaluate(x::MinimumAtom)
  return Base.minimum(evaluate(x.children[1]))
end

# x >= this if minimum(x) = this
# so, x - this will be in the :NonNeg cone
function conic_form(x::MinimumAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    this = Variable()
    objective, constraints = conic_form(this, unique_constr)
    expr = x.children[1] - this
    _, expr_constraints = conic_form(expr, unique_constr)
    append!(constraints, expr_constraints)
    new_constraint = ConicConstr([objective], :NonNeg, [get_vectorized_size(x)])
    push!(constraints, new_constraint)
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end

minimum(x::AbstractExpr) = MinimumAtom(x)
