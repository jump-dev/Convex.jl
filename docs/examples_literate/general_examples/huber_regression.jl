# # Huber regression

# This example can be found here: <https://web.stanford.edu/~boyd/papers/pdf/cvx_applications.pdf>.
# Here we set `big_example = false` to only generate a small example which takes less time to run.
big_example = false
if big_example
    n = 300
    number_tests = 50
else
    n = 50
    number_tests = 10
end

# Generate data for Huber regression.
using Random
Random.seed!(1);
number_samples = round(Int, 1.5 * n);
beta_true = 5 * randn(n);
X = randn(n, number_samples);
Y = zeros(number_samples);
v = randn(number_samples);

#-

## Generate data for different values of p.
## Solve the resulting problems.
using Convex, SCS, Distributions
lsq_data = zeros(number_tests);
huber_data = zeros(number_tests);
prescient_data = zeros(number_tests);
p_vals = range(0, stop = 0.15, length = number_tests);
for i in 1:length(p_vals)
    p = p_vals[i]
    ## Generate the sign changes.
    factor = 2 * rand(Binomial(1, 1 - p), number_samples) .- 1
    Y = factor .* X' * beta_true + v

    ## Form and solve a standard regression problem.
    beta = Variable(n)
    fit = norm(beta - beta_true) / norm(beta_true)
    cost = norm(X' * beta - Y)
    prob = minimize(cost)
    solve!(prob, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    lsq_data[i] = evaluate(fit)

    ## Form and solve a prescient regression problem,
    ## i.e., where the sign changes are known.
    cost = norm(factor .* (X' * beta) - Y)
    solve!(
        minimize(cost),
        MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0),
    )
    prescient_data[i] = evaluate(fit)

    ## Form and solve the Huber regression problem.
    cost = sum(huber(X' * beta - Y, 1))
    solve!(
        minimize(cost),
        MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0),
    )
    huber_data[i] = evaluate(fit)
end

#-
using Plots

plot(p_vals, huber_data, label = "Huber", xlabel = "p", ylabel = "Fit")
plot!(p_vals, lsq_data, label = "Least squares")
plot!(p_vals, prescient_data, label = "Prescient")
#-

## Plot the relative reconstruction error for Huber and prescient regression,
## zooming in on smaller values of p.
indices = findall(p_vals .<= 0.08);
plot(p_vals[indices], huber_data[indices], label = "Huber")
plot!(p_vals[indices], prescient_data[indices], label = "Prescient")
