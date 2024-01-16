# This atom computes dot(sort(x), sort(w)), where w is constant. For example, if
# w = [1 1 1 0 0 0 ... 0], it computes the sum of the 3 largest elements of x
mutable struct DotSortAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    w::Value

    function DotSortAtom(x::AbstractExpr, w::Value)
        if sign(x) == ComplexSign()
            error(
                "[DotSortAtom] argument should be real instead it is $(sign(x))",
            )
        elseif length(w) != length(x)
            error("[DotSortAtom] x and w must be the same size")
        end
        return new((x,), (1, 1), reshape(w, length(x)))
    end
end

head(io::IO, ::DotSortAtom) = print(io, "dotsort")

function Base.sign(x::DotSortAtom)
    if all(x.w .>= 0)
        return sign(x.children[1])
    elseif all(x.w .<= 0)
        # FIXME: shouldn't this be -sign?
        return sign(x.children[1])
    else
        return NoSign()
    end
end

function monotonicity(x::DotSortAtom)
    if all(x.w .>= 0)
        return (Nondecreasing(),)
    end
    # FIXME: if all(x.w .<= 0) isn't it (Nonincreasing(),)?
    return (NoMonotonicity(),)
end

curvature(::DotSortAtom) = ConvexVexity()

function evaluate(x::DotSortAtom)
    y = only(x.children)
    return sum(a * b for (a, b) in zip(sort(vec(evaluate(y))), sort(x.w)))
end

function new_conic_form!(context::Context{T}, x::DotSortAtom) where {T}
    y = only(x.children)
    if size(y, 1) > 1 && size(y, 2) > 1
        y = vec(y)
    end
    μ, ν, e = Variable(size(y)), Variable(size(y)), ones(T, size(y))
    add_constraint!(context, y * x.w' <= e * ν' + μ * e')
    return conic_form!(context, sum(μ) + sum(ν))
end

dotsort(a::AbstractExpr, b::Value) = DotSortAtom(a, b)
dotsort(b::Value, a::AbstractExpr) = DotSortAtom(a, b)
