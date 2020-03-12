#############################################################################
# expressions.jl
# Defines AbstractExpr, which is subtyped by all atoms
# Each type which subtypes AbstractExpr (Variable and Constant being exceptions)
# must have:
#
## head::Symbol                  -- a symbol such as :norm, :+ etc
## children::(AbstractExpr,)     -- The expressions on which the current expression
##                               -- is operated
## id_hash::UInt64               -- identifier hash, can be a hash of children
##                                  or a unique identifier of the object
## size::(Int, Int)          -- size of the resulting expression
#
# Constants and variables do not have children.
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
## conic_form!: TODO: Fill this in after conic_form! is stable
#
#############################################################################

import Base.sign, Base.size, Base.length, Base.lastindex, Base.ndims, Base.convert, Base.axes

### Abstract types
abstract type AbstractExpr end
abstract type Constraint end

const hashaa_seed = hash(AbstractExpr)

id_hash(x::AbstractExpr) = getfield(x, :id_hash)
head(x::AbstractExpr) = getfield(x, :head)
combinehash(x1, x2) = xor(x1*hashaa_seed, x2)

Base.hash(ex::AbstractExpr, h::UInt) = xor(id_hash(ex), hash(head(ex), h))::UInt

function Base.hash(A::Array{<:AbstractExpr}, h::UInt)
    h = xor(hashaa_seed, hash(size(A), h))
    h ⊻= reduce(combinehash, map(id_hash, A))
    return h
end



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

function length(x::AbstractExpr)
    return prod(x.size)
end

### User-defined Unions
const Value = Union{Number, AbstractArray}
const ValueOrNothing = Union{Value, Nothing}
const AbstractExprOrValue = Union{AbstractExpr, Value}

convert(::Type{AbstractExpr}, x::Value) = Constant(x)
convert(::Type{AbstractExpr}, x::AbstractExpr) = x

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
lastindex(x::AbstractExpr) = length(x)

axes(x::AbstractExpr) = (Base.OneTo(size(x, 1)), Base.OneTo(size(x, 2)))
axes(x::AbstractExpr, n::Integer) = axes(x)[n]
lastindex(x::AbstractExpr, n::Integer) = last(axes(x, n))

@deprecate get_vectorized_size(x::AbstractExpr) length(x)
