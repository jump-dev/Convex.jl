# # Trade-off curves
using Random
Random.seed!(1)
m = 25;
n = 10;
A = randn(m, n);
b = randn(m, 1);

#-

using Convex, SCS, LinearAlgebra


gammas = exp10.(range(-4, stop=2, length=100));

x_values = zeros(n, length(gammas));
x = Variable(n);
for i=1:length(gammas)
    cost = sumsquares(A*x - b) + gammas[i]*norm(x,1);
    problem = minimize(cost, [norm(x, Inf) <= 1]);
    solve!(problem, () -> SCS.Optimizer(verbose=0));
    x_values[:,i] = evaluate(x);
end

#-

# Plot the regularization path.

using Plots
plot(title = "Entries of x vs lambda", xaxis=:log, xlabel="lambda", ylabel="x" )
for i = 1:n
    plot!(gammas, x_values[i,:], label="x$i")
end
plot!()
