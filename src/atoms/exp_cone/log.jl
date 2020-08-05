#############################################################################
# log.jl
# natural logarithm of an logression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.log

### Logarithm

mutable struct LogAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function LogAtom(x::AbstractExpr)
        if sign(x)==ComplexSign()
            error("The argument should be real but it's instead complex")
        else
            children = (x,)
            return new(:log, hash(children), children, x.size)
        end
    end
end

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

function template(e::LogAtom, context::Context)
     # log(z) \geq x  <=>    (x,ones(),z) \in ExpCone
     z = e.children[1]
     y = Constant(ones(size(z)))
     x = Variable(size(z))
     objective = template(x, context)
     add_constraints_to_context(ExpConstraint(x, y, z), context)
    return objective
end
