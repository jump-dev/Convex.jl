# deprecate these soon
norm_inf(x::AbstractExpr) = maximum(abs(x))

norm_1(x::AbstractExpr) = sum(abs(x))

norm_fro(x::AbstractExpr) = LinearAlgebra.norm2(vec(x))

"""
    norm(x::AbstractExpr, p::Real=2)

Computes the `p`-norm `‖x‖ₚ = (∑ᵢ |xᵢ|^p)^(1/p)` of a vector expression `x`.

Matrices are vectorized (i.e., `norm(x)` is the same as `norm(vec(x))`.)

This function uses specialized methods for `p=1, 2, Inf`. For `p > 1` otherwise,
this function uses the procedure documented at
[`rational_to_socp.pdf`](https://github.com/jump-dev/Convex.jl/raw/master/docs/supplementary/rational_to_socp.pdf),
based on the paper "Second-order cone programming" by F. Alizadeh and D. Goldfarb,
Mathematical Programming, Series B, 95:3-51, 2001.
"""
function LinearAlgebra.norm(x::AbstractExpr, p::Real = 2)
    if size(x, 2) > 1
        x = vec(x)
    end
    # x is a vector
    if p == 1
        return norm_1(x)
    elseif p == 2
        return LinearAlgebra.norm2(x)
    elseif p == Inf
        return norm_inf(x)
    elseif p > 1
        # TODO: allow tolerance in the rationalize step
        return rationalnorm(x, rationalize(Int, float(p)))
    else
        error("vector p-norms not defined for p < 1")
    end
end
