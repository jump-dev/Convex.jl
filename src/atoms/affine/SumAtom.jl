mutable struct SumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    SumAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::SumAtom) = print(io, "sum")

Base.sign(x::SumAtom) = sign(x.children[1])

monotonicity(::SumAtom) = (Nondecreasing(),)

curvature(::SumAtom) = ConstVexity()

evaluate(x::SumAtom) = sum(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, A::SumAtom) where {T}
    subobj = conic_form!(context, only(AbstractTrees.children(A)))
    return operate(sum, T, sign(A), subobj)
end

# Dispatch to an internal helper function that handles the dimension argument in
# the same manner as Base, with dims=: denoting a regular sum
Base.sum(x::AbstractExpr; dims = :) = _sum(x, dims)

_sum(x::AbstractExpr, ::Colon) = SumAtom(x)

function _sum(x::AbstractExpr, dimension::Integer)
    if dimension == 1
        return Constant(ones(1, x.size[1]), Positive()) * x
    elseif dimension == 2
        return x * Constant(ones(x.size[2], 1), Positive())
    else
        error("[SumAtom] sum not implemented for `dims=$dimension`")
    end
end
