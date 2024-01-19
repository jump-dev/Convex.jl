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

@deprecate norm_inf(x::AbstractExpr) norm(x, Inf)

@deprecate norm_1(x::AbstractExpr) norm(x, 1)

@deprecate norm_fro(x::AbstractExpr) norm(x, 2)

@deprecate get_vectorized_size(x::AbstractExpr) length(x)

@deprecate operatornorm(x::AbstractExpr) LinearAlgebra.opnorm(x)

function Base.in(x::AbstractExpr, y::Symbol)
    if !(y in (:semidefinite, :SDP))
        error("Set $y not understood")
    end
    @warn(
        "Using `x in $y` to construct a semidefinite constraint is " *
        "deprecated. Use `isposdef(x)` or `x ⪰ 0` instead.",
        maxlog = 1,
    )
    return isposdef(x)
end
