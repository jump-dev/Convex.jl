using Convex, LinearAlgebra, SCS
# https://github.com/baggepinnen/EuclideanDistanceMatrices.jl/blob/4a2a1662846ea949786aabab323a3cdb21255dfa/src/EuclideanDistanceMatrices.jl#L57
function complete_distmat(D, W, λ=sqrt(count(W .== 0)))
    @assert all(==(1), diag(W)) "The diagonal is always observed and equal to 0. Make sure the diagonal of W is true"
    @assert all(iszero, diag(D)) "The diagonal of D is always 0"
    n = size(D, 1)
    x = -1/(n + sqrt(n))
    y = -1/sqrt(n)
    V = [fill(y, 1, n-1); fill(x, n-1,n-1) + I(n-1)]
    e = ones(n)
    G = Convex.Variable((n-1, n-1))
    B = V*G*V'
    E = diag(B)*e' + e*diag(B)' - 2*B
    problem = Convex.maximize(tr(G)- λ * norm(vec(W .* (E - D))), [G ∈ :SDP])
    Convex.solve2!(problem, () -> SCS.Optimizer(verbose=false))
    Int(problem.status) == 1 || @error problem.status
    B  = Convex.evaluate(B)
    D2 = diag(B)*e' + e*diag(B)' - 2*B
    @info "Data fidelity (norm(W .* (D-D̃))/norm(D))", (norm(W .* (D-D2))/norm(D))
    s  = svd(B)
    D2, s
end
using Distances, Serialization
# P = randn(2,50)
# D = pairwise(SqEuclidean(), P)
# W = rand(size(D)...) .> 0.3 # Create a random mask
# W = (W + W') .> 0           # It makes sense for the mask to be symmetric
# W[diagind(W)] .= true
# D0 = W .* D               # Remove missing entries

# serialize("D0W.jls", (D0, W))

(D0, W) = deserialize(joinpath(@__DIR__, "D0W.jls"))

@info "With USE_SPARSE3"
Convex.USE_SPARSE() = true
Convex.USE_SPARSE2() = true
Convex.USE_SPARSE3() = true
complete_distmat(D0, W);
@time D2, S = complete_distmat(D0, W);
# @profview complete_distmat(D0, W);

# error()
@info "With USE_SPARSE2"
Convex.USE_SPARSE() = true
Convex.USE_SPARSE2() = true
Convex.USE_SPARSE3() = false

complete_distmat(D0, W);
@time D2, S = complete_distmat(D0, W);
# @profview complete_distmat(D0, W);

@info "With USE_SPARSE"
Convex.USE_SPARSE() = true
Convex.USE_SPARSE2() = false
Convex.USE_SPARSE3() = false
complete_distmat(D0, W);
@time D2, S = complete_distmat(D0, W);
nothing
# @show (norm(D-D2)/norm(D))
# @show (norm(W .* (D-D2))/norm(D))

# @profview complete_distmat(D0, W);
