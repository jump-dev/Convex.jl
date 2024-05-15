# # Robust approximate fitting
# Section 6.4.2
# Boyd & Vandenberghe "Convex Optimization"
# Original by Lieven Vandenberghe
# Adapted for Convex by Joelle Skaf - 10/03/05
#
# Adapted for Convex.jl by Karanveer Mohan and David Zeng - 26/05/14
# Original CVX code and plots here:
# <http://web.cvxr.com/cvx/examples/cvxbook/Ch06_approx_fitting/html/fig6_15.html>
#
# Consider the least-squares problem:
#       minimize $\|(A + tB)x - b\|_2$
# where $t$ is an uncertain parameter in [-1,1]
# Three approximate solutions are found:
#
#   1. nominal optimal (that is, letting t=0)
#   2. stochastic robust approximation:
#           minimize $\mathbb{E}\|(A+tB)x - b\|_2$
#      assuming $u$ is uniformly distributed on [-1,1].
#      (reduces to minimizing $\mathbb{E} \|(A+tB)x-b\|^2 = \|A*x-b\|^2  + x^TPx$
#        where $P = \mathbb{E}(t^2) B^TB = (1/3) B^TB$ )
#   3. worst-case robust approximation:
#           minimize $\mathrm{sup}_{-1\leq u\leq 1} \|(A+tB)x - b\|_2$
#      (reduces to minimizing $\max\{\|(A-B)x - b\|_2, \|(A+B)x - b\|_2\}$ ).
#
using Convex, LinearAlgebra, SCS

# Input Data
m = 20;
n = 10;
A = randn(m, n);
(U, S, V) = svd(A);
S = diagm(exp10.(range(-1, stop = 1, length = n)));
A = U[:, 1:n] * S * V';

B = randn(m, n);
B = B / norm(B);

b = randn(m, 1);
x = Variable(n)

# ## Case 1: nominal optimal solution
p = minimize(norm(A * x - b, 2))
solve!(p, SCS.Optimizer; silent_solver = true)
#-
x_nom = evaluate(x)

# ## Case 2: stochastic robust approximation
P = 1 / 3 * B' * B;
p = minimize(square(pos(norm(A * x - b))) + quadform(x, Symmetric(P)))
solve!(p, SCS.Optimizer; silent_solver = true)
#-
x_stoch = evaluate(x)

# ## Case 3: worst-case robust approximation
p = minimize(max(norm((A - B) * x - b), norm((A + B) * x - b)))
solve!(p, SCS.Optimizer; silent_solver = true)

#-
x_wc = evaluate(x)

# ## Plots
# Here we plot the residuals.
parvals = range(-2, stop = 2, length = 100);

errvals(x) = [norm((A + parvals[k] * B) * x - b) for k in eachindex(parvals)]
errvals_ls = errvals(x_nom)
errvals_stoch = errvals(x_stoch)
errvals_wc = errvals(x_wc)

using Plots
plot(parvals, errvals_ls, label = "Nominal problem")
plot!(parvals, errvals_stoch, label = "Stochastic Robust Approximation")
plot!(parvals, errvals_wc, label = "Worst-Case Robust Approximation")
plot!(
    title = "Residual r(u) vs a parameter u for three approximate solutions",
    xlabel = "u",
    ylabel = "r(u) = ||A(u)x-b||_2",
)
