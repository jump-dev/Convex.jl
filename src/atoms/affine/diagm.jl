#############################################################################
# diagm.jl
# Converts a vector of size n into an n x n diagonal
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import LinearAlgebra.diagm, LinearAlgebra.Diagonal

mutable struct DiagMatrixAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function DiagMatrixAtom(x::AbstractExpr)
        (num_rows, num_cols) = x.size

        if num_rows == 1
            sz = num_cols
        elseif num_cols == 1
            sz = num_rows
        else
            throw(
                ArgumentError(
                    "Only vectors are allowed for diagm/Diagonal. Did you mean to use diag?",
                ),
            )
        end

        children = (x,)
        return new(children, (sz, sz))
    end
end

head(io::IO, ::DiagMatrixAtom) = print(io, "diagm")

function sign(x::DiagMatrixAtom)
    return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::DiagMatrixAtom)
    return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::DiagMatrixAtom)
    return ConstVexity()
end

function evaluate(x::DiagMatrixAtom)
    return Diagonal(vec(evaluate(x.children[1])))
end

function diagm((d, x)::Pair{<:Integer,<:AbstractExpr})
    d == 0 || throw(ArgumentError("only the main diagonal is supported"))
    return DiagMatrixAtom(x)
end
Diagonal(x::AbstractExpr) = DiagMatrixAtom(x)
diagm(x::AbstractExpr) = DiagMatrixAtom(x)

function _conic_form!(context::Context{T}, x::DiagMatrixAtom) where {T}
    obj = conic_form!(context, only(children(x)))

    sz = x.size[1]
    I = collect(1:sz+1:sz*sz)
    J = collect(1:sz)
    V = one(T)
    coeff = GBMatrix{T, T}(I, J, V, sz * sz, sz)
    # coeff = sparse(, 1:sz, one(T),

    return operate(add_operation, T, sign(x), coeff, obj)
end
