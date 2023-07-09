import Base.Broadcast.broadcasted

struct QolElemAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QolElemAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size != y.size
            error(
                "elementwise quad over lin must take two arguments of the same size",
            )
        end
        children = (x, y)
        return new(:qol_elem, hash(children), children, x.size)
    end
end

function sign(q::QolElemAtom)
    return Positive()
end

function monotonicity(q::QolElemAtom)
    return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

function curvature(q::QolElemAtom)
    return ConvexVexity()
end

function evaluate(q::QolElemAtom)
    return (evaluate(q.children[1]) .^ 2) ./ evaluate(q.children[2])
end

function conic_form!(context::Context{T}, q::QolElemAtom) where {T}
    sz = q.children[1].size
    t = Variable(sz[1], sz[2])
    t_obj = conic_form!(context, t)
    x, y = q.children

    for i in 1:length(q.children[1])
        add_constraint!(
            context,
            SOCConstraint(y[i] + t[i], y[i] - t[i], 2 * x[i]),
        )
    end
    add_constraint!(context, y >= 0)
    return t_obj
end

qol_elementwise(x::AbstractExpr, y::AbstractExpr) = QolElemAtom(x, y)

function broadcasted(::typeof(^), x::AbstractExpr, k::Int)
    return k == 2 ? QolElemAtom(x, constant(ones(x.size[1], x.size[2]))) :
           error("raising variables to powers other than 2 is not implemented")
end

invpos(x::AbstractExpr) = QolElemAtom(constant(ones(x.size[1], x.size[2])), x)
function broadcasted(::typeof(/), x::Value, y::AbstractExpr)
    return dotmultiply(constant(x), invpos(y))
end
function /(x::Value, y::AbstractExpr)
    return size(y) == (1, 1) ? MultiplyAtom(constant(x), invpos(y)) :
           error("cannot divide by a variable of size $(size(y))")
end
sumsquares(x::AbstractExpr) = square(norm2(x))

function square(x::AbstractExpr)
    if sign(x) == ComplexSign()
        error(
            "Square of complex number is not DCP. Did you mean square_modulus?",
        )
    else
        QolElemAtom(x, constant(ones(x.size[1], x.size[2])))
    end
end
