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

model = problem.model

function get_info(model)
    n_constraints = 0
    n_elements = 0
    n_variables = MOI.get(model, MOI.NumberOfVariables())
    n_constants = 0
    for (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        n_constraints += MOI.get(problem.model, MOI.NumberOfConstraints{F,S}())

        for inds in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            set = MOI.get(model, MOI.ConstraintSet(), inds)
            n_elements += MOI.dimension(set)
            f = MOI.get(model, MOI.ConstraintFunction(), inds)
            n_constants += count_size(f)
        end
    end
    return (; n_variables, n_constraints, n_elements, n_constants)
end

function count_size(f::MOI.VectorAffineFunction)
    return length(f.terms) + length(f.constants)
end

function count_size(x)
    @show typeof(x)
    @show propertynames(x)
    return 0
end
