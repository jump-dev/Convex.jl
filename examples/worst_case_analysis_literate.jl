## Generate data for worst-case risk analysis.
srand(2);
n = 5;
r = abs(randn(n, 1))/15;
Sigma = 0.9 * rand(n, n) - 0.15;
Sigma_nom = Sigma' * Sigma;
Sigma_nom -= (maximum(Sigma_nom) - 0.9)

#-

## Form and solve portfolio optimization problem.
## Here we minimize risk while requiring a 0.1 return.
using Convex, SCS
set_default_solver(SCSSolver(verbose=0));
w = Variable(n);
ret = dot(r, w);
risk = sum(quad_form(w, Sigma_nom));
problem = minimize(risk, [sum(w) == 1, ret >= 0.1, norm(w, 1) <= 2])
solve!(problem);
wval = vec(evaluate(w))

#-

## Form and solve worst-case risk analysis problem.
Sigma = Semidefinite(n);
Delta = Variable(n, n);
risk = sum(quad_form(wval, Sigma));
problem = maximize(risk, [Sigma == Sigma_nom + Delta,
                    diag(Delta) == 0,
                    abs(Delta) <= 0.2,
                    Delta == Delta']);
solve!(problem);
println("standard deviation = ", round(sqrt(wval' * Sigma_nom * wval), 2));
println("worst-case standard deviation = ", round(sqrt(evaluate(risk)), 2));
println("worst-case Delta = ");
println(round(evaluate(Delta), 2));

