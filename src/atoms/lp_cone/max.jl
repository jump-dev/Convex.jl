#############################################################################
# max.jl
# Return the maximum of the two arguments. Operates elementwise over arrays.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import Base.max

# TODO: This can easily be extended to work
### Max Atom
struct MaxAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function MaxAtom(x::AbstractExpr, y::AbstractExpr)
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
        return new(:max, hash(children), children, sz)
    end
end

function sign(x::MaxAtom)
    sign_one = sign(x.children[1])
    sign_two = sign(x.children[2])
    if sign_one == Positive() || sign_two == Positive()
        return Positive()
    elseif sign_one == Negative() && sign_two == Negative()
        return Negative()
    else
        return sign_one + sign_two
    end
end

# The monotonicity
function monotonicity(x::MaxAtom)
    return (Nondecreasing(), Nondecreasing())
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::MaxAtom)
    return ConvexVexity()
end

function evaluate(x::MaxAtom)
    return max.(evaluate(x.children[1]), evaluate(x.children[2]))
end


function template(x::MaxAtom, context::Context)
    t = Variable(x.size)
    for child in children(x)
        add_constraints_to_context(t >= child, context)
    end
    return template(t, context)

end

max(x::AbstractExpr, y::AbstractExpr) = MaxAtom(x, y)
max(x::AbstractExpr, y::Value) = max(x, constant(y))
max(x::Value, y::AbstractExpr) = max(constant(x), y)
pos(x::AbstractExpr) = max(x, constant(0, Positive()))
hinge_loss(x::AbstractExpr) = pos(1 - x)
