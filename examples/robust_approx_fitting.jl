# Section 6.4.2
# Boyd & Vandenberghe "Convex Optimization"
# Original by Lieven Vandenberghe
# Adapted for CVX by Joelle Skaf - 10/03/05
# Adapted for CVX.jl by Karanveer Mohan and David Zeng - 26/05/14
#
# Consider the least-squares problem:
#       minimize ||(A + tB)x - b||_2
# where t is an uncertain parameter in [-1,1]
# Three approximate solutions are found:
#   1- nominal optimal (i.e. letting t=0)
#   2- stochastic robust approximation:
#           minimize E||(A+tB)x - b||_2
#      assuming u is uniformly distributed on [-1,1] )
#      (reduces to minimizing E ||(A+tB)x-b||^2 = ||A*x-b||^2  + x^TPx
#        where P = E(t^2) B^TB = (1/3) B^TB )
#   3- worst-case robust approximation:
#           minimize sup{-1<=u<=1} ||(A+tB)x - b||_2)
#      (reduces to minimizing max{||(A-B)x - b||_2, ||(A+B)x - b||_2} )

# Input Data
m = 20;
n = 10;
A = randn(m,n);
(U,S,V) = svd(A);
S = diagm(logspace(-1, 1, n));
A = U[:, 1:n] * S * V';

B = randn(m, n);
B = B / Base.norm(B);

b = randn(m, 1);
x = Variable(n)

# Case 1: Nominal optimal solution
p = minimize(norm(A * x - b, 2))
solve!(p)
x_nom = x.evaluate()

# Case 2: Stochastic robust approximation
P = 1 / 3 * B' * B;
p = minimize(square_pos(norm(A * x - b)) + quad_form(x, P))
solve!(p)
x_stoch = x.evaluate()

# Case 3: Worst-case robust approximation
p = minimize(max(norm((A - B) * x - b), norm((A + B) * x - b)))
solve!(p)
x_wc = x.evaluate()

# plot residuals
novals = 100;
parvals = linspace(-2, 2, novals);

errvals_ls = [];
errvals_stoch = [];
errvals_wc = [];

for k=1:novals
  errvals_ls = [errvals_ls, Base.norm((A + parvals[k] * B) * x_nom - b)];
  errvals_stoch = [errvals_stoch, Base.norm((A + parvals[k] * B) * x_stoch - b)];
  errvals_wc = [errvals_wc, Base.norm((A + parvals[k] * B) * x_wc - b)];
end

# Plots. You'll need Gaston and gnuplot installed. Otherwise,
# you might have to change the code a bit to work with other plotting
# libraries
using Gaston
Gaston.set_terminal("x11")
plot(parvals, errvals_ls, "color","blue", parvals, errvals_stoch, "color", "black", parvals, errvals_wc, "color","red")
