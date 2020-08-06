import Base.Broadcast.broadcasted

struct QolElemAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function QolElemAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size != y.size
            error("elementwise quad over lin must take two arguments of the same size")
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
    return (evaluate(q.children[1]).^2) ./ evaluate(q.children[2])
end

function conic_form!(q::QolElemAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, q)
        sz = q.children[1].size
        t = Variable(sz[1], sz[2])
        qol_objective = conic_form!(t, unique_conic_forms)
        x, y = q.children
        conic_form!(SOCElemConstraint(y + t, y - t, 2 * x), unique_conic_forms)
        # add implicit constraint y >= 0
        conic_form!(y >= 0, unique_conic_forms)
        cache_conic_form!(unique_conic_forms, q, qol_objective)
    end
    return get_conic_form(unique_conic_forms, q)
end

qol_elementwise(x::AbstractExpr, y::AbstractExpr) = QolElemAtom(x, y)

broadcasted(::typeof(^),x::AbstractExpr,k::Int) = k==2 ? QolElemAtom(x, Constant(ones(x.size[1], x.size[2]))) : error("raising variables to powers other than 2 is not implemented")

invpos(x::AbstractExpr) = QolElemAtom(Constant(ones(x.size[1], x.size[2])), x)
broadcasted(::typeof(/), x::Value, y::AbstractExpr) = DotMultiplyAtom(constant(x), invpos(y))
/(x::Value, y::AbstractExpr) = size(y) == (1,1) ? MultiplyAtom(constant(x), invpos(y)) : error("cannot divide by a variable of size $(size(y))")
sumsquares(x::AbstractExpr) = square(norm2(x))

function square(x::AbstractExpr)
    if sign(x) == ComplexSign()
        error("Square of complex number is not DCP. Did you mean square_modulus?")
    else
        QolElemAtom(x, Constant(ones(x.size[1], x.size[2])))
    end
end
