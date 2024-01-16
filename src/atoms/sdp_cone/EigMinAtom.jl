mutable struct EigMinAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EigMinAtom(x::AbstractExpr)
        if size(x, 1) != size(x, 2)
            error("[EigMinAtom] eigmin can only be applied to a square matrix.")
        end
        return new((x,), (1, 1))
    end
end

head(io::IO, ::EigMinAtom) = print(io, "eigmin")

Base.sign(::EigMinAtom) = NoSign()

monotonicity(::EigMinAtom) = (Nondecreasing(),)

curvature(::EigMinAtom) = ConcaveVexity()

function evaluate(x::EigMinAtom)
    return LinearAlgebra.eigmin(evaluate(x.children[1]))
end

LinearAlgebra.eigmin(x::AbstractExpr) = EigMinAtom(x)

function new_conic_form!(context::Context, x::EigMinAtom)
    A = only(x.children)
    m, n = size(A)
    t = Variable()
    add_constraint!(context, A - t * LinearAlgebra.I(n) ⪰ 0)
    return conic_form!(context, t)
end
