import LinearAlgebra.norm
export norm_inf, norm, norm_1

# deprecate these soon
norm_inf(x::AbstractExpr) = maximum(abs(x))
norm_1(x::AbstractExpr) = sum(abs(x))
norm_fro(x::AbstractExpr) = norm2(vec(x))

# behavior of norm should be consistent with julia:
# * vector norms for vectors
# * operator norms for matrices
function norm(x::AbstractExpr, p::Real=2)
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
        # TODO: After the deprecation period, allow this again but make it consistent with
        # LinearAlgebra, i.e. make norm(x, p) for x a matrix the same as norm(vec(x), p).
        Base.depwarn("""
                     `norm(x, p)` for matrices will in the future be equivalent to
                     `norm(vec(x), p)`. Use `opnorm(x, p)` for the Julia 0.6 behavior of
                     computing the operator norm for matrices.
                     """, :norm)
        return opnorm(x, p)
    end
end

if isdefined(LinearAlgebra, :vecnorm) # deprecated but defined
    import LinearAlgebra: vecnorm
end
Base.@deprecate vecnorm(x::AbstractExpr, p::Real=2) norm(vec(x), p)
