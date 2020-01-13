# This file adds a new syntax to call elementwise operations on Convex expressions
# because overloading broadcast no longer works in Julia v0.6

import LinearAlgebra.dot

# multiplication
dot(::typeof(*)) = applyDotMultiplyAtom
applyDotMultiplyAtom(x, y) = broadcast(*, x, y)

# division
dot(::typeof(/)) = applyDotDivideAtom
applyDotDivideAtom(x, y) = broadcast(/, x, y)

# exponentiation
dot(::typeof(^)) = applyDotExpAtom
applyDotExpAtom(x, y) = broadcast(^, x, y)

# addition and subtraction are no problem since + and .+ should do the same thing

# alternative syntax

# dot(::typeof(*), x::Variable, y::Variable) = broadcast(*, x, y)
# dot(::typeof(/), x::Variable, y::Variable) = broadcast(/, x, y)
# dot(::typeof(+), x::Variable, y::Variable) = broadcast(+, x, y)
# dot(::typeof(-), x::Variable, y::Variable) = broadcast(-, x, y)
# dot(::typeof(^), x::Variable, y::Variable) = broadcast(^, x, y)
