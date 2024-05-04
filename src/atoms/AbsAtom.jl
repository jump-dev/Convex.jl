# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct AbsAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    AbsAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::AbsAtom) = print(io, "abs")

Base.sign(::AbsAtom) = Positive()

monotonicity(x::AbsAtom) = (Nondecreasing() * sign(x.children[1]),)

curvature(::AbsAtom) = ConvexVexity()

evaluate(x::AbsAtom) = abs.(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, A::AbsAtom) where {T}
    x = only(A.children)
    t = Variable(size(x))
    t_obj = conic_form!(context, t)
    if iscomplex(x)
        tape = conic_form!(context, x)
        re_tape = real(tape)
        im_tape = imag(tape)
        if re_tape isa SPARSE_VECTOR
            @assert im_tape isa SPARSE_VECTOR
            return abs.(re_tape + im * im_tape)
        end
        re_safs = MOI.Utilities.scalarize(to_vaf(re_tape))
        im_safs = MOI.Utilities.scalarize(to_vaf(im_tape))
        for (re, im, ti) in zip(re_safs, im_safs, t_obj.variables)
            f = MOI.Utilities.operate(vcat, T, ti, re, im)
            MOI.add_constraint(context.model, f, MOI.SecondOrderCone(3))
        end
    else
        add_constraint!(context, t >= x)
        add_constraint!(context, t >= -x)
    end
    return t_obj
end

Base.abs(x::AbstractExpr) = AbsAtom(x)

Base.abs2(x::AbstractExpr) = square(abs(x))
