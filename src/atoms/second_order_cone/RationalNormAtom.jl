# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

# Handles the k-norms for k > 1 where k is rational. Reduces the
# k-norm constraint to at most 2d ceil(log2(n + m)) + O(d) second
# order cone constraints. Here d is the dimension of the problem and k
# is converted to a rational with value k = n / m. The procedure for
# this is based on the paper "Second-order cone programming" by
# F. Alizadeh and D. Goldfarb, Mathematical Programming, Series B,
# 95:3-51, 2001, and is documented in the pdf available at
# https://github.com/jump-dev/Convex.jl/raw/master/docs/supplementary/rational_to_socp.pdf

mutable struct RationalNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    k::Rational{Int64}

    function RationalNormAtom(x::AbstractExpr, k::Rational{Int})
        if k < 1
            error("[RationalNormAtom] p-norms not defined for p < 1. Got $k")
        end
        return new((x,), (1, 1), k)
    end
end

head(io::IO, ::RationalNormAtom) = print(io, "rationalnorm")

Base.sign(::RationalNormAtom) = Positive()

function monotonicity(x::RationalNormAtom)
    return (sign(x.children[1]) * Nondecreasing(),)
end

curvature(::RationalNormAtom) = ConvexVexity()

function evaluate(x::RationalNormAtom)
    return sum(abs.(evaluate(x.children[1])) .^ x.k)^(1 / x.k)
end

function new_conic_form!(context::Context{T}, x::RationalNormAtom) where {T}
    v = conic_form!(context, only(x.children))
    d = length(only(x.children)) + 1
    t = Variable()
    t_obj = conic_form!(context, t)
    f = operate(vcat, T, sign(x), t_obj, v)
    MOI_add_constraint(context.model, f, MOI.NormCone(Float64(x.k), d))
    return t_obj
end

function rationalnorm(x::AbstractExpr, k::Rational{Int})
    if sign(x) == ComplexSign()
        row, col = size(x)
        if !(row == 1 || col == 1)
            error("[RationalNormAtom] not defined for complex matrices")
        end
        return RationalNormAtom(abs(x), k)
    end
    return RationalNormAtom(x, k)
end
