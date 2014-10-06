#############################################################################
# maximum.jl
# Compute the maximum value of an array.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.maximum
export maximum

### Maximum Atom
type MaximumAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function MaximumAtom(x::AbstractExpr)
    children = (x,)
    return new(:maximum, hash(children), children, (1, 1))
  end
end

function sign(x::MaximumAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::MaximumAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MaximumAtom)
  return ConvexVexity()
end

function evaluate(x::MaximumAtom)
  return Base.maximum(evaluate(x.children[1]))
end

# x <= this if maximum(x) = this
# so, this - x will be in the :NonNeg cone
function conic_form(x::MaximumAtom, unique_constr)
  if !((x.head, x.children_hash) in unique_constr)
    this = Variable()
    objective, constraints = conic_form(this, unique_constr)
    expr = this - x.children[1]
    _, expr_constraints = conic_form(expr, unique_constr)
    append!(constraints, expr_constraints)
    new_constraint = ConicConstr([objective], :NonNeg, [get_vectorized_size(x)])
    push!(constraints, new_constraint)
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end

maximum(x::AbstractExpr) = MaximumAtom(x)
