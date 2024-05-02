# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct DiagMatrixAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function DiagMatrixAtom(x::AbstractExpr)
        num_rows, num_cols = x.size
        if num_rows == 1
            return new((x,), (num_cols, num_cols))
        elseif num_cols == 1
            return new((x,), (num_rows, num_rows))
        else
            msg = "[DiagMatrixAtom] only vectors are allowed for `diagm(x)` and `Diagonal(x). Did you mean to use `diag(x, 0)`?"
            throw(ArgumentError(msg))
        end
    end
end

head(io::IO, ::DiagMatrixAtom) = print(io, "diagm")

Base.sign(x::DiagMatrixAtom) = sign(x.children[1])

monotonicity(::DiagMatrixAtom) = (Nondecreasing(),)

curvature(::DiagMatrixAtom) = ConstVexity()

function evaluate(x::DiagMatrixAtom)
    return LinearAlgebra.Diagonal(vec(evaluate(x.children[1])))
end

function new_conic_form!(context::Context{T}, x::DiagMatrixAtom) where {T}
    I = 1:(x.size[1]+1):x.size[1]^2
    coeff = create_sparse(T, I, 1:x.size[1], one(T), x.size[1]^2, x.size[1])
    obj = conic_form!(context, x.children[1])
    return operate(add_operation, T, sign(x), coeff, obj)
end
