function set_with_size(::Type{MOI.Nonnegatives}, sz::Tuple{Int,Int})
    return MOI.Nonnegatives(prod(sz))
end

head(io::IO, ::MOI.Nonnegatives) = print(io, "≥")

function is_feasible(f, ::MOI.Nonnegatives, tol)
    return all(f .>= -tol)
end

function vexity(vex, ::MOI.Nonnegatives)
    if vex  == ConvexVexity()
        return NotDcp()
    elseif vex == ConcaveVexity()
        return ConvexVexity()
    end
    return vex
end

function _promote_size(lhs::AbstractExpr, rhs::AbstractExpr)
    if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
        error(
            "Cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))",
        )
    end
    if lhs.size == rhs.size || lhs.size == (1, 1)
        sz = rhs.size
        if lhs.size == (1, 1) && rhs.size != (1, 1)
            lhs = lhs * ones(rhs.size)
        end
    elseif rhs.size == (1, 1)
        sz = lhs.size
        if rhs.size == (1, 1) && lhs.size != (1, 1)
            rhs = rhs * ones(lhs.size)
        end
    else
        error(
            "Cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)",
        )
    end
    return lhs, rhs
end

function Base.:>=(lhs::AbstractExpr, rhs::AbstractExpr)
    lhs, rhs = _promote_size(lhs, rhs)
    return GenericConstraint{MOI.Nonnegatives}(lhs - rhs)
end

Base.:>=(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))

Base.:>=(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)
