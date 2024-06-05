mutable struct RootDetAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    RootDetAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::RootDetAtom) = print(io, "rootdet")

Base.sign(::RootDetAtom) = NoSign()

monotonicity(::RootDetAtom) = (NoMonotonicity(),)

curvature(::RootDetAtom) = ConcaveVexity()

function evaluate(x::RootDetAtom)
    n = size(only(AbstractTrees.children(x)), 1)
    return LinearAlgebra.det(evaluate(x.children[1]))^(1 / n)
end

function new_conic_form!(context::Context{T}, x::RootDetAtom) where {T}
    t = conic_form!(context, Variable())
    A = only(AbstractTrees.children(x))
    f = operate(vcat, T, sign(x), t, conic_form!(context, A))
    MOI_add_constraint(context.model, f, MOI.RootDetConeSquare(size(A, 1)))
    return t
end
