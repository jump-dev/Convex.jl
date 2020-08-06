#############################################################################
# min.jl
# Return the minimum of the two arguments. Operates elementwise over arrays.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.min

# TODO: This can easily be extended to work
### Min Atom
struct MinAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function MinAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error("Both the arguments should be real instead they are $(sign(x)) and $(sign(y))")
        else
            if x.size == y.size
                sz = x.size
            elseif x.size == (1, 1)
                sz = y.size
            elseif y.size == (1, 1)
                sz = x.size
            else
                error("Got different sizes for x as $(x.size) and y as $(y.size)")
            end
        end

        children = (x, y)
        return new(:min, hash(children), children, sz)
    end
end

function sign(x::MinAtom)
    sign_one = sign(x.children[1])
    sign_two = sign(x.children[2])
    if sign_one == Negative() || sign_two == Negative()
        return Negative()
    elseif sign_one == Positive() && sign_two == Positive()
        return Positive()
    else
        return sign_one + sign_two
    end
end

# The monotonicity
function monotonicity(x::MinAtom)
    return (Nondecreasing(), Nondecreasing())
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MinAtom)
    return ConcaveVexity()
end

function evaluate(x::MinAtom)
    return min.(evaluate(x.children[1]), evaluate(x.children[2]))
end


function conic_form!(x::MinAtom, context::UniqueConicForms)
    t = Variable(x.size[1], x.size[2])
    for child in x.children
        add_constraints_to_context(t <= child, context)
    end
    return template(t, context)
end


min(x::AbstractExpr, y::AbstractExpr) = MinAtom(x, y)
min(x::AbstractExpr, y::Value) = min(x, constant(y))
min(x::Value, y::AbstractExpr) = min(constant(x), y)
neg(x::AbstractExpr) = max(-x, Constant(0, Positive()))
