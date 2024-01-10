mutable struct QolElemAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QolElemAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size != y.size
            error(
                "elementwise quad over lin must take two arguments of the same size",
            )
        end
        return new((x, y), x.size)
    end
end

head(io::IO, ::QolElemAtom) = print(io, "qol_elem")

Base.sign(::QolElemAtom) = Positive()

function monotonicity(q::QolElemAtom)
    return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

curvature(::QolElemAtom) = ConvexVexity()

function evaluate(q::QolElemAtom)
    return (evaluate(q.children[1]) .^ 2) ./ evaluate(q.children[2])
end

function new_conic_form!(context::Context{T}, q::QolElemAtom) where {T}
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

function Base.Broadcast.broadcasted(::typeof(^), x::AbstractExpr, k::Int)
    if k != 2
        error("raising variables to powers other than 2 is not implemented")
    end
    return QolElemAtom(x, constant(ones(x.size[1], x.size[2])))
end

function Base.Broadcast.broadcasted(
    ::typeof(Base.literal_pow),
    ::typeof(^),
    x::AbstractExpr,
    ::Val{k},
) where {k}
    return Base.Broadcast.broadcasted(^, x, k)
end

invpos(x::AbstractExpr) = QolElemAtom(constant(ones(x.size[1], x.size[2])), x)

function Base.Broadcast.broadcasted(::typeof(/), x::Value, y::AbstractExpr)
    return dotmultiply(constant(x), invpos(y))
end

function Base.:/(x::Value, y::AbstractExpr)
    if size(y) != (1, 1)
        error("cannot divide by a variable of size $(size(y))")
    end
    return MultiplyAtom(constant(x), invpos(y))
end

sumsquares(x::AbstractExpr) = square(norm2(x))

function square(x::AbstractExpr)
    if sign(x) == ComplexSign()
        error(
            "Square of complex number is not DCP. Did you mean square_modulus?",
        )
    end
    return QolElemAtom(x, constant(ones(x.size[1], x.size[2])))
end
