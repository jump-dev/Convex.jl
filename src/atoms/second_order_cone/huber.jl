export huber

struct HuberAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    M::Real

    function HuberAtom(x::AbstractExpr, M::Real)
        if sign(x) == ComplexSign()
            error("argument must be real")
        elseif M <= 0
            error("huber parameter must by a positive scalar")
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

function conic_form!(x::HuberAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        c = x.children[1]
        s = Variable(c.size)
        n = Variable(c.size)

        # objective given by s.^2 + 2 * M * |n|
        objective = conic_form!(square(s) + 2 * x.M * abs(n), unique_conic_forms)
        conic_form!(c == s + n, unique_conic_forms)

        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

huber(x::AbstractExpr, M::Real=1.0) = HuberAtom(x, M)
