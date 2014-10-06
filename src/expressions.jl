#############################################################################
# expressions.jl
# Defines AbstractExpr, which is subtyped by all atoms
# Each type which subtypes AbstractExpr (Variable and Constant being exceptions)
# must have:
#
## head::Symbol                  -- a symbol such as :vecnorm, :+ etc
## children::(AbstractExpr,)     -- The expressions on which the current expression
##                               -- is operated
## children_hash::Uint64         -- hash of all the children
## size::(Int64, Int64)          -- size of the resulting expression
#
# To eliminate as many edge cases as possibles, Constant and Variables don't have
# children, or children_hash as fields. They do have an id, value, size and vexity
# as fields.
#
# In addition, each atom must implement the following functions:
## sign: returns the sign of the result of the expression
## monotonicity: The monotonicity of the arguments with respect to the function
##      i.e if the argument is nondecreasing, will the function be nonincreasing
##      or nondecreasing? eg. negate(x) will have Nonincreasing monotonicity
## evaluate: Evaluates the value of the expression, assuming the problem has been
##           solved.
## curvature: If h(x)=f∘g(x), then (for single variable calculus)
##      h''(x) = g'(x)^T f''(g(x)) g'(x) + f'(g(x))g''(x)
##      curvature refers to the curvature of the first term.
##      We then use this curvature to find vexity of h (see vexity function below)
## conic_form: TODO: Fill this in after conic_form is stable
#
#############################################################################

import Base.sign, Base.size, Base.endof, Base.ndims
export AbstractExpr
export vexity, sign, size, evaluate, monotonicity, curvature
export conic_form
export endof, ndims
export Value, ValueOrNothing
export get_vectorized_size

### Abstract type
abstract AbstractExpr

# If h(x)=f∘g(x), then (for single variable calculus)
# h''(x) = g'(x)^T f''(g(x)) g'(x) + f'(g(x))g''(x)
# We calculate the vexity according to this
function vexity(x::AbstractExpr)
  monotonicities = monotonicity(x)
  vex = curvature(x)
  for i = 1:length(x.children)
    vex += monotonicities[i] * vexity(x.children[i])
  end
  return vex
end

# This function should never be reached
function monotonicity(x::AbstractExpr)
  error("monotonicity not implemented for $(x.head).")
end

# This function should never be reached
function curvature(x::AbstractExpr)
  error("curvature not implemented for $(x.head).")
end

# This function should never be reached
function evaluate(x::AbstractExpr)
  error("evaluate not implemented for $(x.head).")
end

# This function should never be reached
function sign(x::AbstractExpr)
  error("sign not implemented for $(x.head).")
end

function size(x::AbstractExpr)
  return x.size
end

### User-defined Unions
Value = Union(Number, AbstractArray)
ValueOrNothing = Union(Value, Nothing)

### Indexing Utilities
endof(x::AbstractExpr) = x.size[1] * x.size[2]

function size(x::AbstractExpr, dim::Integer)
  if dim < 1
    error("dimension out of range")
  elseif dim > 2
    return 1
  else
    return size(x)[dim]
  end
end

ndims(x::AbstractExpr) = 2

get_vectorized_size(x::AbstractExpr) = reduce(*, size(x))
