mutable struct QuadOverLinAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QuadOverLinAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size[2] != 1 && y.size != (1, 1)
            error("quad over lin arguments must be a vector and a scalar")
        end
        return new((x, y), (1, 1))
    end
end

head(io::IO, ::QuadOverLinAtom) = print(io, "qol")

Base.sign(::QuadOverLinAtom) = Positive()

function monotonicity(q::QuadOverLinAtom)
    return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

curvature(::QuadOverLinAtom) = ConvexVexity()

function evaluate(q::QuadOverLinAtom)
    x = evaluate(q.children[1])
    return x' * x / evaluate(q.children[2])
end

function new_conic_form!(context::Context, q::QuadOverLinAtom)
    t = Variable()
    x, y = q.children
    add_constraint!(context, SOCConstraint(y + t, y - t, 2 * x))
    add_constraint!(context, y >= 0)
    return conic_form!(context, t)
end

quadoverlin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)
