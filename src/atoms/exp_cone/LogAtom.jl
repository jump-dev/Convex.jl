mutable struct LogAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[LogAtom] the argument should be real but it's instead complex",
            )
        end
        return new((x,), x.size)
    end
end

head(io::IO, ::LogAtom) = print(io, "log")

Base.sign(::LogAtom) = NoSign()

monotonicity(::LogAtom) = (Nondecreasing(),)

curvature(::LogAtom) = ConcaveVexity()

evaluate(x::LogAtom) = log.(evaluate(x.children[1]))

Base.log(x::AbstractExpr) = LogAtom(x)

function new_conic_form!(context::Context, e::LogAtom)
    # log(z) \geq x  <=> (x,ones(),z) \in ExpCone
    z = e.children[1]
    y = Constant(ones(size(z)))
    x = Variable(size(z))
    add_constraint!(context, ExponentialConeConstraint(vcat(x, y, z)))
    return conic_form!(context, x)
end
