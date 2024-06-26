# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct ImaginaryAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    ImaginaryAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::ImaginaryAtom) = print(io, "imag")

Base.sign(::ImaginaryAtom) = NoSign()

monotonicity(::ImaginaryAtom) = (Nondecreasing(),)

curvature(::ImaginaryAtom) = ConstVexity()

evaluate(x::ImaginaryAtom) = imag.(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, x::ImaginaryAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(imag, T, sign(x), obj)
end
