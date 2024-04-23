mutable struct ExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[ExpAtom] the argument should be real but it's instead complex",
            )
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
    # exp(x) \leq z  <=>  (x,1,z) \in ExpCone
    x = e.children[1]
    m, n = size(x)
    z = Variable(m, n)
    for i in 1:m, j in 1:n
        f = vcat(x[i, j], 1, z[i, j])
        add_constraint!(context, GenericConstraint{MOI.ExponentialCone}(f))
    end
    return conic_form!(context, z)
end
