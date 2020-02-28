# # SVM with L^1 regularization
## Generate data for SVM classifier with L1 regularization.
using Random
Random.seed!(3);
n = 20;
m = 1000;
TEST = m;
DENSITY = 0.2;
beta_true = randn(n,1);
idxs = randperm(n)[1:round(Int, (1-DENSITY)*n)];
beta_true[idxs] .= 0
offset = 0;
sigma = 45;
X = 5 * randn(m, n);
Y = sign.(X * beta_true .+ offset .+ sigma * randn(m,1));
X_test = 5 * randn(TEST, n);

#-

## Form SVM with L1 regularization problem.
using Convex, SCS, ECOS

beta = Variable(n);
v = Variable();
loss = sum(pos(1 - Y .* (X*beta - v)));
reg = norm(beta, 1);

## Compute a trade-off curve and record train and test error.
TRIALS = 100
train_error = zeros(TRIALS);
test_error = zeros(TRIALS);
lambda_vals = exp10.(range(-2, stop=0, length=TRIALS);)
beta_vals = zeros(length(beta), TRIALS);
for i = 1:TRIALS
    lambda = lambda_vals[i];
    problem = minimize(loss/m + lambda*reg);
    solve!(problem, () -> ECOS.Optimizer(verbose=0));
    ## solve!(problem, SCS.Optimizer(verbose=0,linear_solver=SCS.Direct, eps=1e-3))
    train_error[i] = sum(float(sign.(X*beta_true .+ offset) .!= sign.(evaluate(X*beta - v))))/m;
    test_error[i] = sum(float(sign.(X_test*beta_true .+ offset) .!= sign.(evaluate(X_test*beta - v))))/TEST;
    beta_vals[:, i] =  evaluate(beta);
end

#-

# Plot the train and test error over the trade-off curve.
using Plots
plot(lambda_vals, train_error, label="Train error");
plot!(lambda_vals, test_error, label="Test error");
plot!(xscale=:log, yscale=:log, ylabel="errors", xlabel="lambda")

#-

# Plot the regularization path for beta.

plot()
for i = 1:n
    plot!(lambda_vals, vec(beta_vals[i,:]), label="beta$i")
end
plot!(xscale=:log, ylabel="betas", xlabel="lambda")
