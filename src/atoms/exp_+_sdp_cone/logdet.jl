import LinearAlgebra.logdet

struct LogDetAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function LogDetAtom(x::AbstractExpr)
        children = (x,)
        return new(:logdet, hash(children), children, (1, 1))
    end
end

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


function template(x::LogDetAtom, context::Context{T}) where {T}
    # the object we want the logdet of. Should be a PSD matrix, but may not be a `AbstractVariable` itself.
    A = only(children(x)) 

    # We vectorize and take the upper triangle
    v = vec_triu(A)

    # ensure symmetry
    add_constraints_to_context(v == vec_tril(A), context)

    # We pass to MOI
    X = template(v, context)

    t = template(Variable(), context)
    f = operate(vcat, T, t, [1], X)
    side_dimension = size(only(children(x)), 1)

    set = MOI.LogDetConeTriangle(side_dimension)

    MOI_add_constraint(context.model, f,set)
    return t
end

logdet(x::AbstractExpr) = LogDetAtom(x)
