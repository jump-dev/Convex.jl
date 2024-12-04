# # Entropy Maximization
#
# Here is a constrained entropy maximization problem:
#
# ```math
# \begin{array}{ll}
#     \text{maximize}   & -\sum_{i=1}^n x_i \log x_i \\
#     \text{subject to} & \mathbf{1}' x = 1 \\
#                   & Ax \leq b
# \end{array}
# ```
#
# where $x \in \mathbf{R}^n$ is our optimization variable and $A \in \mathbf{R}^{m \times n}, b \in \mathbf{R}^{m}$.
#
# To solve this, we can simply use the `entropy` operation Convex.jl provides.

using Convex, SCS

n = 25;
m = 15;
A = randn(m, n);
b = rand(m, 1);

x = Variable(n);
problem = maximize(entropy(x), sum(x) == 1, A * x <= b)
solve!(problem, SCS.Optimizer; silent = true)

#-

evaluate(x)
