#############################################################################
# exp.jl
# e raised to the power of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.exp

### Exponential

struct ExpAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function ExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("The argument should be real but it's instead complex")
        else
            children = (x,)
            return new(:exp, hash(children), children, x.size)
        end
    end
end

function sign(x::ExpAtom)
    return Positive()
end

function monotonicity(x::ExpAtom)
    return (Nondecreasing(),)
end

function curvature(x::ExpAtom)
    return ConvexVexity()
end

function evaluate(x::ExpAtom)
    return exp.(evaluate(x.children[1]))
end

exp(x::AbstractExpr) = ExpAtom(x)


function template(e::ExpAtom, context::Context{T}) where {T}
    # exp(x) \leq z  <=>  (x,ones(),z) \in ExpCone
    x = e.children[1]
    y = Constant(ones(size(x)))
    z = Variable(size(x))
    add_constraints_to_context(ExpConstraint(x, y, z), context)
    return template(z, context)
end
