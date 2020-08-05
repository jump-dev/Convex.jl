using ECOS
using Convex
using MathOptInterface
const MOI = MathOptInterface
using JuMP
using Test

# Generate fake data matrix
function gen_data(m, n, k)
    return (10 * rand(m, k) *  2 *rand(k, n))
end

function gen_masks(A, holdout)
    training_mask = rand(size(A)...) .> 1 - holdout
    validation_mask = .!training_mask
    return training_mask, validation_mask
end

function alternating_minimization(f, A, M, Y_init, k, MAX_ITERS)
    m, n = size(A)

    X = Variable(m, k)
    Y = Variable(k, n)

    objective = (norm(vec(M .* A) - vec(M .* (X*Y)), 2)
                 + γ1 * norm(vec(X), 2)
                 + γ2 * norm(vec(Y), 1))

    constraints =  [X*Y >= ϵ]

    problem = minimize(objective, constraints)

    Y.value = Y_init
    for i in 1:MAX_ITERS
        fix!(Y)
        res = f(problem)
        free!(Y)
        fix!(X)
        res = f(problem)
        free!(X)
    end

    return problem, X, Y
end

function JuMP_setup!(model, X, Y, A, M, k)
    m, n = size(A)

    @variable(model, t1 >= 0)
    @variable(model, t2 >= 0)
    @variable(model, Y_abs[1:n*k] >= 0)
    @constraint(model, vcat(t1, vec(M .* A - M .* (X*Y))) ∈ SecondOrderCone())
    @constraint(model, vcat(t2, vec(X)) ∈ SecondOrderCone())
    @constraint(model, vec(Y) .<= Y_abs )
    @constraint(model, -vec(Y) .<= Y_abs )
    @objective(model, Min, t1 + γ1 * t2 + γ2*sum(Y_abs))
    @constraint(model, X*Y .>= ϵ)
end

function alternating_minimization_JuMP(A, M, Y_init, k, MAX_ITERS)
    m, n = size(A)

    function Y_fix_model(Y)
        model = Model(() -> ECOS.Optimizer(verbose=false))
        @variable(model, X[1:m, j=1:k])
        JuMP_setup!(model, X, Y, A, M, k)
        return X, model
    end

    function X_fix_model(X)
        model = Model(() -> ECOS.Optimizer(verbose=false))
        @variable(model, Y[1:k, j=1:n])
        JuMP_setup!(model, X, Y, A, M, k)
        return Y, model
    end

    local Xval, model
    Yval = Y_init
    for i in 1:MAX_ITERS
        X, model = Y_fix_model(Yval)
        JuMP.optimize!(model)
        Xval = value.(X)

        Y, model = X_fix_model(Xval)
        JuMP.optimize!(model)
        Yval = value.(Y)
    end

    return model, Xval, Yval
end

const γ1 = 1.0
const γ2 = 1.0
const ϵ = 0.0001
MAX_ITERS = 2

m, n, k = 150, 150, 5
holdout = 0.80

A = gen_data(m, n, k)
Mt, Mv = gen_masks(A, holdout)

Y_init = rand(k,n)
# @info "Running with classic `Convex.solve!`..." (MAX_ITERS,m,n,k)
# @time p1, X1, Y1 = alternating_minimization(A, Mt, Y_init, k, MAX_ITERS) do problem
#     solve!(problem, () -> ECOS.Optimizer(verbose=false))
# end

Convex.USE_SPARSE() = true
# recompile
alternating_minimization(A, Mt, Y_init, k, MAX_ITERS) do problem
    solve2!(problem, ECOS.Optimizer(verbose=false))
end

@info "Running with `Convex.solve2!`..." (MAX_ITERS,m,n,k, Convex.USE_SPARSE())
@time p2c, X2c, Y2c = alternating_minimization(A, Mt, Y_init, k, MAX_ITERS) do problem
    solve2!(problem, ECOS.Optimizer(verbose=false))
end

# Convex.USE_SPARSE() = false
# # recompile
# alternating_minimization(A, Mt, Y_init, k, MAX_ITERS) do problem
#     solve2!(problem, ECOS.Optimizer(verbose=false))
# end

# @info "Running with `Convex.solve2!`..." (MAX_ITERS,m,n,k, Convex.USE_SPARSE())
# @time p2, X2, Y2 = alternating_minimization(A, Mt, Y_init, k, MAX_ITERS) do problem
#     solve2!(problem, ECOS.Optimizer(verbose=false))
# end

@info "Running with JuMP..." (MAX_ITERS,m,n,k)
@time model, X3, Y3 = alternating_minimization_JuMP(A, Mt, Y_init, k, MAX_ITERS);

@testset "Same results" begin
    # @test evaluate(X1) ≈ evaluate(X2) atol=1e-2 rtol=1e-2
    # @test evaluate(Y1) ≈ evaluate(Y2) atol=1e-2 rtol=1e-2

    # @test evaluate(X1) ≈ X3 atol=1e-2 rtol=1e-2
    # @test evaluate(Y1) ≈ Y3 atol=1e-2 rtol=1e-2

    # @test evaluate(X2c) ≈ evaluate(X2) atol=1e-2 rtol=1e-2
    # @test evaluate(Y2c) ≈ evaluate(Y2) atol=1e-2 rtol=1e-2

    @test evaluate(X2c) ≈ X3 atol=1e-2 rtol=1e-2
    @test evaluate(Y2c) ≈ Y3 atol=1e-2 rtol=1e-2
end
