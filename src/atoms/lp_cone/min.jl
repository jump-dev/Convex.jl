mutable struct MinAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MinAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error(
                "Both the arguments should be real instead they are $(sign(x)) and $(sign(y))",
            )
        end
        sz = if x.size == y.size
            x.size
        elseif x.size == (1, 1)
            y.size
        elseif y.size == (1, 1)
            x.size
        else
            error("Got different sizes for x as $(x.size) and y as $(y.size)")
        end
        return new((x, y), sz)
    end
end

head(io::IO, ::MinAtom) = print(io, "min")

function Base.sign(x::MinAtom)
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

monotonicity(::MinAtom) = (Nondecreasing(), Nondecreasing())

curvature(::MinAtom) = ConcaveVexity()

evaluate(x::MinAtom) = min.(evaluate(x.children[1]), evaluate(x.children[2]))

function new_conic_form!(context::Context, x::MinAtom)
    t = Variable(x.size[1], x.size[2])
    p = maximize(t, (t <= child for child in x.children)...)
    return conic_form!(context, p)
end

Base.min(x::AbstractExpr, y::AbstractExpr) = MinAtom(x, y)

Base.min(x::AbstractExpr, y::Value) = min(x, constant(y))

Base.min(x::Value, y::AbstractExpr) = min(constant(x), y)

neg(x::AbstractExpr) = max(-x, Constant(0, Positive()))
