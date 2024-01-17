mutable struct MaxAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MaxAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error(
                "[MaxAtom] both the arguments should be real instead they are $(sign(x)) and $(sign(y))",
            )
        end
        sz = if x.size == y.size
            x.size
        elseif x.size == (1, 1)
            y.size
        elseif y.size == (1, 1)
            x.size
        else
            error(
                "[MaxAtom] got different sizes for x as $(x.size) and y as $(y.size)",
            )
        end
        return new((x, y), sz)
    end
end

head(io::IO, ::MaxAtom) = print(io, "max")

function Base.sign(x::MaxAtom)
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

monotonicity(::MaxAtom) = (Nondecreasing(), Nondecreasing())

curvature(::MaxAtom) = ConvexVexity()

evaluate(x::MaxAtom) = max.(evaluate(x.children[1]), evaluate(x.children[2]))

function new_conic_form!(context::Context, x::MaxAtom)
    t = Variable(x.size)
    t_obj = conic_form!(context, t)
    for child in x.children
        add_constraint!(context, t >= child)
    end
    return t_obj
end

Base.max(x::AbstractExpr, y::AbstractExpr) = MaxAtom(x, y)

Base.max(x::AbstractExpr, y::Value) = max(x, constant(y))

Base.max(x::Value, y::AbstractExpr) = max(constant(x), y)

pos(x::AbstractExpr) = max(x, Constant(0, Positive()))

hinge_loss(x::AbstractExpr) = pos(1 - x)
