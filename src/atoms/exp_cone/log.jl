#############################################################################
# log.jl
# natural logarithm of an logression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.log

### Logarithm

mutable struct LogAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("The argument should be real but it's instead complex")
        else
            children = (x,)
            return new(children, x.size)
        end
    end
end
head(io::IO, ::LogAtom) = print(io, "log")

function sign(x::LogAtom)
    return NoSign()
end

function monotonicity(x::LogAtom)
    return (Nondecreasing(),)
end

function curvature(x::LogAtom)
    return ConcaveVexity()
end

function evaluate(x::LogAtom)
    return log.(evaluate(x.children[1]))
end

log(x::AbstractExpr) = LogAtom(x)

function conic_form!(context::Context, e::LogAtom)
    # log(z) \geq x  <=>    (x,ones(),z) \in ExpCone
    z = e.children[1]
    y = Constant(ones(size(z)))
    x = Variable(size(z))
    objective = conic_form!(context, x)
    add_constraint!(context, ExpConstraint(x, y, z))
    return objective
end
