struct HuberAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    M::Real

    function HuberAtom(x::AbstractExpr, M::Real)
        if sign(x) == ComplexSign()
            error("Arguemt must be real")
        elseif M <= 0
            error("Huber parameter must by a positive scalar")
        end
        children = (x,)
        return new(:huber, hash((children, M)), children, x.size, M)
    end
end

function sign(x::HuberAtom)
    return Positive()
end

function monotonicity(x::HuberAtom)
    return (Nondecreasing() * sign(x.children[1]),)
end

function curvature(x::HuberAtom)
    return ConvexVexity()
end

function evaluate(x::HuberAtom)
    c = evaluate(x.children[1])
    for i in 1:length(c)
        if c[i] <= x.M
            c[i] = c[i]^2
        else
            c[i] = 2*x.M*c[i] - x.M^2
        end
    end
    return c
end

function template(x::HuberAtom, context::Context)
    c = x.children[1]
    s = Variable(c.size)
    n = Variable(c.size)

    # objective given by s.^2 + 2 * M * |n|
    objective = template(square(s) + 2 * x.M * abs(n), context)
    add_constraints_to_context(c == s + n, context)
    return objective
end

huber(x::AbstractExpr, M::Real=1.0) = HuberAtom(x, M)
