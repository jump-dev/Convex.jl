#############################################################################
# abs.jl
# Absolute value of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.abs, Base.abs2

### Absolute Value

struct AbsAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function AbsAtom(x::AbstractExpr)
        children = (x,)
        return new(:abs, hash(children), children, x.size)
    end
end

function sign(x::AbsAtom)
    return Positive()
end

function monotonicity(x::AbsAtom)
    return (Nondecreasing() * sign(x.children[1]),)
end

function curvature(x::AbsAtom)
    return ConvexVexity()
end

function evaluate(x::AbsAtom)
    return abs.(evaluate(x.children[1]))
end

function conic_form!(context::Context, A::AbsAtom)
    x = only(A.children)

    t = Variable(size(x))
    t_obj = conic_form!(context, t)

    if sign(x) == ComplexSign()
        for i in 1:length(vec(t))
            add_constraints_to_context(
                t[i] >= norm2([real(x[i]); imag(x[i])]),
                context,
            )
        end
    else
        add_constraints_to_context(t >= x, context)
        add_constraints_to_context(t >= -x, context)
    end
    return t_obj
end

abs(x::AbstractExpr) = AbsAtom(x)
abs2(x::AbstractExpr) = square(abs(x))

## Alternate version: no atom, just a DCPPromise:

# # A lightweight atom
# struct DCPPromiseWrapper <: AbstractExpr
#     x::AbstractVariable
#     sign::Sign
#     monotonicity::Monotonicity
#     curvature::Vexity
#     evaluate::Any
# end

# function DCPPromiseWrapper(x::AbstractVariable;
#         sign::Sign = sign(x),
#         monotonicity::Monotonicity = monotonicity(x),
#         curvature::Vexity = curvature(x),
#         evaluate = () -> evaluate(x))
#         DCPPromise(tuple(x), sign, monotonicity, curvature, evaluate)
# end

# sign(d::DCPPromiseWrapper) = d.sign
# monotonicity(d::DCPPromiseWrapper) = d.monotonicity
# curvature(d::DCPPromiseWrapper) = d.curvature

# # No children! it is a leaf type
# AbstractTrees.children(::DCPPromiseWrapper) = ()

# Base.size(d::DCPPromiseWrapper) = size(d.x)
# evaluate(d::DCPPromiseWrapper) = d.evaluate()

# conic_form!(context::Context, D::DCPPromiseWrapper) = conic_form!(context, D.x)

# function abs(x::AbstractExpr)
#     t = Variable(size(x), Positive())
#     add_constraint!(t, t >= x)
#     add_constraint!(t, t >= -x)
#     return DCPPromise(
#         t,
#         monotonicity = Nondecreasing() * sign(x),
#         curvature = ConvexVexity(),
#         evaluate = () -> abs.(evaluate(x)),
#     )
# end

# Base.abs2(x::AbstractExpr) = square(abs(x))
