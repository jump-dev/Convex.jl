#############################################################################
# exp.jl
# e raised to the power of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.exp

### Exponential

mutable struct ExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error("The argument should be real but it's instead complex")
        else
            children = (x,)
            return new(children, x.size)
        end
    end
end

head(io::IO, ::ExpAtom) = print(io, "exp")

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

function _conic_form!(context::Context{T}, e::ExpAtom) where {T}
    # exp(x) \leq z  <=>  (x,ones(),z) \in ExpCone
    x = e.children[1]
    y = Constant(ones(size(x)))
    z = Variable(size(x))
    add_constraint!(context, ExpConstraint(x, y, z))
    return conic_form!(context, z)
end
