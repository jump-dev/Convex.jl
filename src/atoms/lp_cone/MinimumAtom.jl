mutable struct MinimumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function MinimumAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[MinimumAtom] argument should be real instead it is $(sign(x))",
            )
        end
        return new((x,), (1, 1))
    end
end

head(io::IO, ::MinimumAtom) = print(io, "minimum")

Base.sign(x::MinimumAtom) = sign(x.children[1])

monotonicity(::MinimumAtom) = (Nondecreasing(),)

curvature(::MinimumAtom) = ConcaveVexity()

evaluate(x::MinimumAtom) = minimum(evaluate(x.children[1]))

function new_conic_form!(context::Context, x::MinimumAtom)
    t = Variable()
    add_constraint!(context, t <= x.children[1])
    return conic_form!(context, t)
end

Base.minimum(x::AbstractExpr) = MinimumAtom(x)
