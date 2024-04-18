mutable struct QolElemAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function QolElemAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size != y.size
            error(
                "[QolElemAtom] elementwise quad over lin must take two arguments of the same size",
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
    x, y = q.children
    t = Variable(x.size)
    for i in 1:length(x)
        f = vcat(y[i] + t[i], y[i] - t[i], 2 * x[i])
        add_constraint!(context, GenericConstraint{MOI.SecondOrderCone}(f))
    end
    add_constraint!(context, y >= 0)
    return conic_form!(context, t)
end

qol_elementwise(x::AbstractExpr, y::AbstractExpr) = QolElemAtom(x, y)

function square(x::AbstractExpr)
    if sign(x) == ComplexSign()
        error(
            "square of a complex number is not DCP. Did you mean square_modulus?",
        )
    end
    return QolElemAtom(x, constant(ones(x.size)))
end

sumsquares(x::AbstractExpr) = square(norm2(x))

invpos(x::AbstractExpr) = QolElemAtom(constant(ones(x.size)), x)

function Base.Broadcast.broadcasted(::typeof(/), x::Value, y::AbstractExpr)
    return _dot_multiply(constant(x), invpos(y))
end

function Base.:/(x::Value, y::AbstractExpr)
    if size(y) != (1, 1)
        error("cannot divide by a variable of size $(size(y))")
    end
    return MultiplyAtom(constant(x), invpos(y))
end

function Base.Broadcast.broadcasted(::typeof(^), x::AbstractExpr, k::Int)
    if k != 2
        error("raising variables to powers other than 2 is not implemented")
    end
    return square(x)
end

function Base.Broadcast.broadcasted(
    ::typeof(Base.literal_pow),
    ::typeof(^),
    x::AbstractExpr,
    ::Val{k},
) where {k}
    return Base.Broadcast.broadcasted(^, x, k)
end
