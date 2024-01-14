# TODO(odow): Remove code in a future release

function Base.:<(x::AbstractExpr, y::AbstractExpr)
    @warn(
        "Strict inequality operator `<` has never enforced the strict " *
        "inequality, and instead adds the non-strict inequality `<=`. This " *
        "operator will be removed in the next release.",
        maxlog = 1,
    )
    return x <= y
end

Base.:<(lhs::AbstractExpr, rhs::Value) = <(lhs, constant(rhs))

Base.:<(lhs::Value, rhs::AbstractExpr) = <(constant(lhs), rhs)

function Base.:>(x::AbstractExpr, y::AbstractExpr)
    @warn(
        "Strict inequality operator `>` has never enforced the strict " *
        "inequality, and instead adds the non-strict inequality `>=`. This " *
        "operator will be removed in the next release.",
        maxlog = 1,
    )
    return x >= y
end

Base.:>(lhs::AbstractExpr, rhs::Value) = >(lhs, constant(rhs))

Base.:>(lhs::Value, rhs::AbstractExpr) = >(constant(lhs), rhs)
