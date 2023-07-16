import LinearAlgebra.logdet

mutable struct LogDetAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogDetAtom(x::AbstractExpr)
        children = (x,)
        return new(children, (1, 1))
    end
end

head(io::IO, ::LogDetAtom) = print(io, "logdet")

function sign(x::LogDetAtom)
    return NoSign()
end

function monotonicity(x::LogDetAtom)
    return (NoMonotonicity(),)
end

function curvature(x::LogDetAtom)
    return ConcaveVexity()
end

function evaluate(x::LogDetAtom)
    return log(det(evaluate(x.children[1])))
end

function _conic_form!(context::Context{T}, x::LogDetAtom) where {T}
    # the object we want the logdet of. Should be a PSD matrix, but may not be a `AbstractVariable` itself.
    A = only(children(x))

    # We vectorize and take the upper triangle
    v = vec_triu(A)

    # ensure symmetry
    add_constraint!(context, v == vec_tril(A))

    # We pass to MOI
    X = conic_form!(context, v)

    t = conic_form!(context, Variable())
    f = operate(vcat, T, sign(x), t, SPARSE_VECTOR{T}([one(T)]), X)
    side_dimension = size(only(children(x)), 1)

    set = MOI.LogDetConeTriangle(side_dimension)

    MOI_add_constraint(context.model, f, set)
    return t
end

logdet(x::AbstractExpr) = LogDetAtom(x)
