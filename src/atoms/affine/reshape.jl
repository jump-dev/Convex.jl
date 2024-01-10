# We process all objects on a vectorized level. The `ReshapeAtom` therefore
# is simply in charge of remembering the new size.
mutable struct ReshapeAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ReshapeAtom(x::AbstractExpr, m::Int, n::Int)
        if m * n != length(x)
            error("Cannot reshape expression of size $(x.size) to ($(m), $(n))")
        end
        return new((x,), (m, n))
    end
end

head(io::IO, ::ReshapeAtom) = print(io, "reshape")

Base.sign(x::ReshapeAtom) = sign(x.children[1])

monotonicity(::ReshapeAtom) = (Nondecreasing(),)

curvature(::ReshapeAtom) = ConstVexity()

function evaluate(x::ReshapeAtom)
    val = evaluate(x.children[1])
    if val isa Number
        return val
    end
    return reshape(val, x.size[1], x.size[2])
end

function new_conic_form!(context::Context, A::ReshapeAtom)
    return conic_form!(context, only(AbstractTrees.children(A)))
end

Base.reshape(x::AbstractExpr, m::Int, n::Int) = ReshapeAtom(x, m, n)

Base.vec(x::AbstractExpr) = reshape(x, length(x), 1)
