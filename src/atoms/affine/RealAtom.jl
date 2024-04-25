# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct RealAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    RealAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::RealAtom) = print(io, "real")

function Base.sign(x::RealAtom)
    if sign(x.children[1]) == ComplexSign()
        return NoSign()
    end
    return sign(x.children[1])
end

monotonicity(::RealAtom) = (Nondecreasing(),)

curvature(::RealAtom) = ConstVexity()

evaluate(x::RealAtom) = real.(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, x::RealAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(real, T, sign(x), obj)
end

Base.real(x::AbstractExpr) = RealAtom(x)

Base.real(x::ComplexConstant) = x.real_constant

Base.real(x::Constant) = x
