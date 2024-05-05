# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct RelativeEntropyAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function RelativeEntropyAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error(
                "[RelativeEntropyAtom] both the arguments should be real but these are instead $(sign(x)) and $(sign(y))",
            )
        end
        return new((x, y), (1, 1))
    end
end

head(io::IO, ::RelativeEntropyAtom) = print(io, "relative_entropy")

Base.sign(::RelativeEntropyAtom) = NoSign()

monotonicity(::RelativeEntropyAtom) = (NoMonotonicity(), NoMonotonicity())

curvature(::RelativeEntropyAtom) = ConvexVexity()

function evaluate(e::RelativeEntropyAtom)
    x = vectorize(evaluate(e.children[1]))
    y = vectorize(evaluate(e.children[2]))
    if any(isnan, y)
        return Inf
    end
    out = x .* log.(x ./ y)
    # fix value when x=0:
    # out will only be NaN if x=0, in which case the correct value is 0
    out[isnan.(out)] .= 0
    return sum(out)
end

function new_conic_form!(context::Context{T}, e::RelativeEntropyAtom) where {T}
    u = Variable()
    x, y = e.children
    # Convex.jl has the connvention:
    #   u >= relative_entropy(x,y) = sum_i( x_i log (x_i/y_i) )
    # But MOI has the reversed convention:
    #   MOI.RelativeEntropyCone has the order (u, y, x)
    f = vcat(u, y, x)
    add_constraint!(context, GenericConstraint{MOI.RelativeEntropyCone}(f))
    return conic_form!(context, u)
end

relative_entropy(x::AbstractExpr, y::AbstractExpr) = RelativeEntropyAtom(x, y)

# y*log(x/y)
log_perspective(x::AbstractExpr, y::AbstractExpr) = -relative_entropy(y, x)
