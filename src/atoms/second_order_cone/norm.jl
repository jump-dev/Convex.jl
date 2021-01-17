import LinearAlgebra.norm

# deprecate these soon
norm_inf(x::AbstractExpr) = maximum(abs(x))
norm_1(x::AbstractExpr) = sum(abs(x))
norm_fro(x::AbstractExpr) = norm2(vec(x))

"""
    norm(x::AbstractExpr, p::Real=2)

Computes the `p`-norm `‖x‖ₚ = (∑ᵢ |xᵢ|^p)^(1/p)` of a vector expression `x`.

This function uses specialized methods for `p=1, 2, Inf`. For `p > 1` otherwise,
this function uses the procedure documented at
[`rational_to_socp.pdf`](https://github.com/jump-dev/Convex.jl/raw/master/docs/supplementary/rational_to_socp.pdf),
based on the paper "Second-order cone programming" by F. Alizadeh and D. Goldfarb,
Mathematical Programming, Series B, 95:3-51, 2001.

!!! warning
    For versions of Convex.jl prior to v0.14.0, `norm` on a matrix expression returned
    the operator norm ([`opnorm`](@ref)), which matches Julia v0.6 behavior. This functionality
    was deprecated since Convex.jl v0.8.0, and has been removed. In the future,
    `norm(x, p)` will return `‖vec(x)‖ₚ`, matching the behavior of [`norm`](@ref)
    for numeric matrices.
"""
function LinearAlgebra.norm(x::AbstractExpr, p::Real=2)
    if length(size(x)) <= 1 || minimum(size(x))==1
        # x is a vector
        if p == 1
            return norm_1(x)
        elseif p == 2
            return norm2(x)
        elseif p == Inf
            return norm_inf(x)
        elseif p > 1
            # TODO: allow tolerance in the rationalize step
            return rationalnorm(x, rationalize(Int, float(p)))
        else
            error("vector p-norms not defined for p < 1")
        end
    else
        error("In Convex.jl v0.13 and below, `norm(x, p)` meant `opnorm(x, p)` (but was deprecated since v0.8.0). In the future, `norm(x,p)` for matrices will be equivalent to `norm(vec(x),p)`. This is currently an error to ensure you update your code!")
    end
end
