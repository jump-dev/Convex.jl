function set_with_size(::Type{MOI.Nonpositives}, sz::Tuple{Int,Int})
    return MOI.Nonpositives(prod(sz))
end

head(io::IO, ::MOI.Nonpositives) = print(io, "≤")

function is_feasible(f, ::MOI.Nonpositives, tol)
    return all(f .<= -tol)
end

function vexity(vex, ::MOI.Nonpositives)
    if vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function Base.:<=(lhs::AbstractExpr, rhs::AbstractExpr)
    lhs, rhs = _promote_size(lhs, rhs)
    return GenericConstraint{MOI.Nonpositives}(lhs - rhs)
end

Base.:<=(lhs::AbstractExpr, rhs::Value) = <=(lhs, constant(rhs))

Base.:<=(lhs::Value, rhs::AbstractExpr) = <=(constant(lhs), rhs)
