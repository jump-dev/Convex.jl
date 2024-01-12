mutable struct ExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("The argument should be real but it's instead complex")
        end
        return new((x,), x.size)
    end
end

head(io::IO, ::ExpAtom) = print(io, "exp")

Base.sign(::ExpAtom) = Positive()

monotonicity(::ExpAtom) = (Nondecreasing(),)

curvature(::ExpAtom) = ConvexVexity()

evaluate(x::ExpAtom) = exp.(evaluate(x.children[1]))

Base.exp(x::AbstractExpr) = ExpAtom(x)

function new_conic_form!(context::Context{T}, e::ExpAtom) where {T}
    # exp(x) \leq z  <=>  (x,ones(),z) \in ExpCone
    x = e.children[1]
    y = Constant(ones(size(x)))
    z = Variable(size(x))
    add_constraint!(context, ExponentialConeConstraint(x, y, z))
    return conic_form!(context, z)
end
