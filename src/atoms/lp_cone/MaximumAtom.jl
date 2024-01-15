mutable struct MaximumAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function MaximumAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[MaximumAtom] argument should be real instead it is $(sign(x))",
            )
        end
        return new((x,), (1, 1))
    end
end

head(io::IO, ::MaximumAtom) = print(io, "maximum")

Base.sign(x::MaximumAtom) = sign(x.children[1])

monotonicity(::MaximumAtom) = (Nondecreasing(),)

curvature(::MaximumAtom) = ConvexVexity()

evaluate(x::MaximumAtom) = Base.maximum(evaluate(x.children[1]))

function new_conic_form!(context::Context, x::MaximumAtom)
    t = Variable()
    add_constraint!(context, t >= x.children[1])
    return conic_form!(context, t)
end

Base.maximum(x::AbstractExpr) = MaximumAtom(x)
