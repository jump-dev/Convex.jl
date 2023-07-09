struct QuadOverLinAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QuadOverLinAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size[2] != 1 && y.size != (1, 1)
            error("quad over lin arguments must be a vector and a scalar")
        end
        children = (x, y)
        return new(:qol, hash(children), children, (1, 1))
    end
end

function sign(q::QuadOverLinAtom)
    return Positive()
end

function monotonicity(q::QuadOverLinAtom)
    return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

function curvature(q::QuadOverLinAtom)
    return ConvexVexity()
end

function evaluate(q::QuadOverLinAtom)
    x = evaluate(q.children[1])
    return x' * x / evaluate(q.children[2])
end

function conic_form!(context::Context, q::QuadOverLinAtom)
    t = Variable()
    x, y = q.children
    add_constraints_to_context(SOCConstraint(y + t, y - t, 2 * x), context)
    add_constraints_to_context(y >= 0, context)
    return conic_form!(context, t)
end

quadoverlin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)
