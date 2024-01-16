mutable struct EigMaxAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EigMaxAtom(x::AbstractExpr)
        if size(x, 1) != size(x, 2)
            error("[EigMaxAtom] eigmax can only be applied to a square matrix.")
        end
        return new((x,), (1, 1))
    end
end

head(io::IO, ::EigMaxAtom) = print(io, "eigmax")

Base.sign(::EigMaxAtom) = NoSign()

monotonicity(::EigMaxAtom) = (Nondecreasing(),)

curvature(::EigMaxAtom) = ConvexVexity()

function evaluate(x::EigMaxAtom)
    return LinearAlgebra.eigmax(evaluate(x.children[1]))
end

LinearAlgebra.eigmax(x::AbstractExpr) = EigMaxAtom(x)

function new_conic_form!(context::Context{T}, x::EigMaxAtom) where {T}
    A = only(x.children)
    m, n = size(A)
    t = Variable()
    add_constraint!(context, t * LinearAlgebra.I(n) - A ⪰ 0)
    return conic_form!(context, t)
end
