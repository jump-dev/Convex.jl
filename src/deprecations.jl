# TODO(odow): Remove code in a future release
function Base.:<(::AbstractExpr, ::AbstractExpr)
    return error("Strict inequality `<` has been removed. Use `<=` instead.")
end
Base.:<(lhs::AbstractExpr, rhs::Value) = <(lhs, constant(rhs))
Base.:<(lhs::Value, rhs::AbstractExpr) = <(constant(lhs), rhs)

function Base.:>(::AbstractExpr, ::AbstractExpr)
    return error("Strict inequality `>` has been removed. Use `>=` instead.")
end
Base.:>(lhs::AbstractExpr, rhs::Value) = >(lhs, constant(rhs))
Base.:>(lhs::Value, rhs::AbstractExpr) = >(constant(lhs), rhs)
