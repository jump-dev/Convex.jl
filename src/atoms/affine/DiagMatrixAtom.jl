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

function LinearAlgebra.diagm((d, x)::Pair{<:Integer,<:AbstractExpr})
    if d != 0
        msg = "[DiagMatrixAtom] only the main diagonal is supported. Got `d=$d`"
        throw(ArgumentError(msg))
    end
    return diagm(x)
end

LinearAlgebra.Diagonal(x::AbstractExpr) = diagm(x)

LinearAlgebra.diagm(x::AbstractExpr) = DiagMatrixAtom(x)

function new_conic_form!(context::Context{T}, x::DiagMatrixAtom) where {T}
    I = 1:(x.size[1]+1):x.size[1]^2
    coeff = create_sparse(T, I, 1:x.size[1], one(T), x.size[1]^2, x.size[1])
    obj = conic_form!(context, x.children[1])
    return operate(add_operation, T, sign(x), coeff, obj)
end
