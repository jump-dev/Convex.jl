mutable struct HuberAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    M::Real

    function HuberAtom(x::AbstractExpr, M::Real)
        if sign(x) == ComplexSign()
            error("Argument must be real")
        elseif M <= 0
            error("Huber parameter must by a positive scalar")
        end
        return new((x,), x.size, M)
    end
end

head(io::IO, ::HuberAtom) = print(io, "huber")

Base.sign(::HuberAtom) = Positive()

monotonicity(x::HuberAtom) = (Nondecreasing() * sign(x.children[1]),)

curvature(::HuberAtom) = ConvexVexity()

function evaluate(x::HuberAtom)
    c = evaluate(x.children[1])
    for i in 1:length(c)
        if c[i] <= x.M
            c[i] = c[i]^2
        else
            c[i] = 2 * x.M * c[i] - x.M^2
        end
    end
    return c
end

function new_conic_form!(context::Context, x::HuberAtom)
    c = x.children[1]
    s = Variable(c.size)
    n = Variable(c.size)
    add_constraint!(context, c == s + n)
    # objective given by s.^2 + 2 * M * |n|
    return conic_form!(context, square(s) + 2 * x.M * abs(n))
end

huber(x::AbstractExpr, M::Real = 1.0) = HuberAtom(x, M)
