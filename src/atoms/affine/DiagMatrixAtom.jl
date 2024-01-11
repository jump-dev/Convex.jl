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
            msg = "Only vectors are allowed for diagm/Diagonal. Did you mean to use diag?"
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
        throw(ArgumentError("only the main diagonal is supported"))
    end
    return DiagMatrixAtom(x)
end

LinearAlgebra.Diagonal(x::AbstractExpr) = DiagMatrixAtom(x)

LinearAlgebra.diagm(x::AbstractExpr) = DiagMatrixAtom(x)

function new_conic_form!(context::Context{T}, x::DiagMatrixAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    sz = x.size[1]
    I = collect(1:sz+1:sz*sz)
    J = collect(1:sz)
    V = one(T)
    coeff = create_sparse(T, I, J, V, sz * sz, sz)
    return operate(add_operation, T, sign(x), coeff, obj)
end
