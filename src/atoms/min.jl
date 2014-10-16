#############################################################################
# min.jl
# Return the minimum of the two arguments. Operates elementwise over arrays.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.min
export min, neg

# TODO: This can easily be extended to work
### Min Atom
type MinAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr, AbstractExpr)
  size::(Int64, Int64)

  function MinAtom(x::AbstractExpr, y::AbstractExpr)
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
    return new(:min, hash(children), children, sz)
  end
end

function sign(x::MinAtom)
  sign_one = sign(x.children[1])
  sign_two = sign(x.children[2])
  if sign_one == Negative() || sign_two == Negative()
    return Negative()
  elseif sign_one == Positive() && sign_two == Positive()
    return Positive()
  else
    return sign_one + sign_two
  end
end

# The monotonicity
function monotonicity(x::MinAtom)
  return (Nonincreasing(), Nonincreasing())
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MinAtom)
  return ConcaveVexity()
end

function evaluate(x::MinAtom)
  return Base.min([evaluate(x)], [evaluate(y)])
end

# x >= this and y >= this if min(x, y) = this
function conic_form(x::MinAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    this = Variable(x.size[1], x.size[2])
    objective = conic_form(this, unique_conic_forms)
    for child in x.children
      conic_form(this <= child, unique_conic_forms)
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

min(x::AbstractExpr, y::AbstractExpr) = MinAtom(x, y)
min(x::AbstractExpr, y::Value) = min(x, Constant(y))
min(x::Value, y::AbstractExpr) = min(Constant(x), y)
neg(x::AbstractExpr) = min(x, Constant(0, Negative()))
