# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

#############################################################################
# expressions.jl
# Defines AbstractExpr, which is subtyped by all atoms
# Each type which subtypes AbstractExpr (Variable and Constant being exceptions)
# must have:
#
## children::(AbstractExpr,)     -- The expressions on which the current expression
##                               -- is operated
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
##
## Optional: `head` to define custom printing
#
#############################################################################

abstract type AbstractExpr end

const Value = Union{Number,AbstractArray}

# We commandeer `==` to create a constraint.
# Therefore we define `isequal` to still have a notion of equality
# (Normally `isequal` falls back to `==`, so we need to provide a method).
# All `AbstractExpr` (Constraints are not AbstractExpr's!) are compared by value, except for AbstractVariables, which are compared by `===` (objectid).
function Base.isequal(x::AbstractExpr, y::AbstractExpr)
    if typeof(x) != typeof(y)
        return false
    end
    for i in 1:fieldcount(typeof(x))
        if !isequal(getfield(x, i), getfield(y, i))
            return false
        end
    end
    return true
end

# Define hash consistently with `isequal`
function Base.hash(x::AbstractExpr, h::UInt)
    h = hash(typeof(x), h)
    for i in 1:fieldcount(typeof(x))
        h = hash(getfield(x, i), h)
    end
    return h
end

# If h(x)=f∘g(x), then (for single variable calculus)
# h''(x) = g'(x)^T f''(g(x)) g'(x) + f'(g(x))g''(x)
# We calculate the vexity according to this
function vexity(x::AbstractExpr)
    monotonicities = monotonicity(x)
    vex = curvature(x)
    for i in 1:length(x.children)
        vex += monotonicities[i] * vexity(x.children[i])
    end
    return vex
end

vexity(::Value) = ConstVexity()

evaluate(x) = output(x) # fallback

Base.size(x::AbstractExpr) = x.size

Base.length(x::AbstractExpr) = prod(x.size)

Base.convert(::Type{AbstractExpr}, x::Value) = constant(x)

Base.convert(::Type{AbstractExpr}, x::AbstractExpr) = x

function Base.size(x::AbstractExpr, dim::Integer)
    if dim < 1
        error("dimension out of range")
    elseif dim > 2
        return 1
    end
    return size(x)[dim]
end

Base.ndims(x::AbstractExpr) = 2

Base.lastindex(x::AbstractExpr) = length(x)

Base.axes(x::AbstractExpr) = (Base.OneTo(size(x, 1)), Base.OneTo(size(x, 2)))

Base.axes(x::AbstractExpr, n::Integer) = axes(x)[n]

Base.lastindex(x::AbstractExpr, n::Integer) = last(axes(x, n))
