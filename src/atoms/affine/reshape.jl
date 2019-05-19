import Base.reshape, Base.vec
export reshape, vec
export sign, curvature, monotonicity, evaluate, conic_form!


struct ReshapeAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function ReshapeAtom(x::AbstractExpr, m::Int, n::Int)
        if m * n != length(x)
            error("cannot reshape expression of size $(x.size) to ($(m), $(n))")
        end
        return new(:reshape, objectid(x), (x,), (m, n))
    end
end

function sign(x::ReshapeAtom)
    return sign(x.children[1])
end

function monotonicity(x::ReshapeAtom)
    return (Nondecreasing(),)
end

function curvature(x::ReshapeAtom)
    return ConstVexity()
end

function evaluate(x::ReshapeAtom)
    return reshape(evaluate(x.children[1]), x.size[1], x.size[2])
end

function conic_form!(x::ReshapeAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    return conic_form!(x.children[1], unique_conic_forms)
end


reshape(x::AbstractExpr, m::Int, n::Int) = ReshapeAtom(x, m, n)
vec(x::AbstractExpr) = reshape(x, length(x), 1)
