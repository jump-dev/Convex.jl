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
    size::Tuple{Int, Int}

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

function conic_form!(x::AbsAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        c = x.children[1]
        t = Variable(size(c))
        objective = conic_form!(t, unique_conic_forms)
        if sign(x.children[1]) == ComplexSign()
            for i in 1:length(vec(t))
                conic_form!(t[i]>=norm2([real(c[i]);imag(c[i])]), unique_conic_forms)
            end
        else
            conic_form!(c <= t, unique_conic_forms)
            conic_form!(c >= -t, unique_conic_forms)
        end
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

function template(A::AbsAtom, context::Context)
    x = only(A.children)
    
    t = Variable(size(x))
    t_obj = template(t, context)

    add_constraints_to_context(t >= x, context)
    add_constraints_to_context(t >= -x, context)
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

# template(D::DCPPromiseWrapper, context::Context) = template(D.x, context)


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
