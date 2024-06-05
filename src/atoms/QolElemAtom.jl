# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct QolElemAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QolElemAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size != y.size
            error(
                "[QolElemAtom] elementwise quad over lin must take two arguments of the same size",
            )
        end
        return new((x, y), x.size)
    end
end

head(io::IO, ::QolElemAtom) = print(io, "qol_elem")

Base.sign(::QolElemAtom) = Positive()

function monotonicity(q::QolElemAtom)
    return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

curvature(::QolElemAtom) = ConvexVexity()

function evaluate(q::QolElemAtom)
    return (evaluate(q.children[1]) .^ 2) ./ evaluate(q.children[2])
end

function new_conic_form!(context::Context{T}, q::QolElemAtom) where {T}
    x, y = q.children
    t = Variable(x.size)
    t_tape = conic_form!(context, t)
    y_tape = conic_form!(context, y)
    x_tape = conic_form!(context, x)
    x_fn = MOI.Utilities.scalarize(to_vaf(x_tape))
    y_fn = MOI.Utilities.scalarize(to_vaf(y_tape))
    for (ti, yi, xi) in zip(t_tape.variables, y_fn, x_fn)
        f = MOI.Utilities.operate(vcat, T, ti, 1 / T(2) * yi, xi)
        MOI.add_constraint(context.model, f, MOI.RotatedSecondOrderCone(3))
    end
    return t_tape
end
