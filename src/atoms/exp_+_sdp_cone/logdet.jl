mutable struct LogDetAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    LogDetAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::LogDetAtom) = print(io, "logdet")

Base.sign(::LogDetAtom) = NoSign()

monotonicity(::LogDetAtom) = (NoMonotonicity(),)

curvature(::LogDetAtom) = ConcaveVexity()

evaluate(x::LogDetAtom) = log(LinearAlgebra.det(evaluate(x.children[1])))

function new_conic_form!(context::Context{T}, x::LogDetAtom) where {T}
    # the object we want the logdet of. Should be a PSD matrix, but may not be a
    # `AbstractVariable` itself.
    A = only(AbstractTrees.children(x))
    v = vec_triu(A)
    add_constraint!(context, v == vec_tril(A))
    t = conic_form!(context, Variable())
    X = conic_form!(context, v)
    f = operate(vcat, T, sign(x), t, SPARSE_VECTOR{T}([one(T)]), X)
    MOI_add_constraint(context.model, f, MOI.LogDetConeTriangle(size(A, 1)))
    return t
end

LinearAlgebra.logdet(x::AbstractExpr) = LogDetAtom(x)
