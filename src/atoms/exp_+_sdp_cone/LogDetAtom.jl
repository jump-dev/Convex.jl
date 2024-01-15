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
    t = conic_form!(context, Variable())
    A = only(AbstractTrees.children(x))
    X = conic_form!(context, A)
    u = SPARSE_VECTOR{T}([one(T)])
    f = operate(vcat, T, sign(x), t, u, X)
    MOI_add_constraint(context.model, f, MOI.LogDetConeSquare(size(A, 1)))
    return t
end

LinearAlgebra.logdet(x::AbstractExpr) = LogDetAtom(x)
