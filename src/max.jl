#############################################################################
# max.jl
# Return the maximum of the two arguments. Operates elementwise over arrays.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.max
export max

# TODO: This can easily be extended to work
### Max Atom
type MaxAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function MaxAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size == y.size
      sz = x.size
    elseif x.size == (1, 1)
      sz = y.size
    elseif y.size == (1, 1)
      sz = x.size
    else
      error("Got different sizes for x as $(x.size) and y as $(y.size)")
    end

    children = (x, y)
    return new(:max, hash(children), children, sz)
  end
end

function sign(x::MaxAtom)
  sign_one = sign(x.children[1])
  sign_two = sign(x.children[2])
  if sign_one == Positive() || sign_two == Positive()
    return Positive()
  elseif sign_one == Negative() && sign_two == Negative()
    return Negative()
  else
    return sign_one + sign_two
  end
end

# The monotonicity
function monotonicity(x::MaxAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MaxAtom)
  return ConvexVexity()
end

function evaluate(x::MaxAtom)
  return Base.max([evaluate(x)], [evaluate(y)])
end

# x <= this and y <= this if max(x, y) = this
function conic_form(x::MaxAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    this = Variable(x.size[1], x.size[2])
    objective, constraints = conic_form(this, unique_constr)
    for child in x.children
      expr = this - child
      _, expr_constraints = conic_form(expr, unique_constr)
      append!(constraints, expr_constraints)
    end
    new_constraint = ConicConstr([objective], :NonNeg, [get_vectorized_size(x)])
    push!(constraints, new_constraint)
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end

max(x::AbstractExpr, y::AbstractExpr) = MaxAtom(x)
