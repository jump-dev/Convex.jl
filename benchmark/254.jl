using Convex, Clarabel, JLD2, Downloads
import MathOptInterface as MOI
file = JLD2.load(
    Downloads.download(
        "https://github.com/cossio/cvx_example_data/raw/master/cvx_example.jld2",
    ),
)

m = 2000

let
    Adj = file["Adj"]
    N = file["N"]
    if m < size(Adj, 1)
        Adj = Adj[1:m, :]
        N = N[1:m, :, :]
    end
    x = Variable(2730)
    S, V, T = size(N)
    lN = log.(N)
    M = vec(sum(N[:, :, 2:T]; dims = (2, 3)))
    H = Adj * x
    problem = maximize(
        -sum(logsumexp(lN[:, v, t] - H) for t in 1:T-1, v in 1:V) - dot(M, H),
        [x ≥ -1e2, x ≤ 1e2],
    )

    @time context = Convex.Context(problem, MOI.Utilities.Model{Float64})
end
