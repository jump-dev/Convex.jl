Quick Tutorial
==============

Consider a constrained least squares problem

```math
\begin{aligned}
\begin{array}{ll}
\text{minimize} & \|Ax - b\|_2^2 \\
\text{subject to} & x \geq 0
\end{array}
\end{aligned}
```

with variable $x\in \mathbf{R}^{n}$, and problem data
$A \in \mathbf{R}^{m \times n}$, $b \in \mathbf{R}^{m}$.

This problem can be solved in Convex.jl as follows: :

```@repl
# Make the Convex.jl module available
using Convex, SCS

# Generate random problem data
m = 4;  n = 5
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
problem = minimize(sumsquares(A * x - b), [x >= 0])

# Solve the problem by calling solve!
solve!(problem, () -> SCS.Optimizer(verbose=false))

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval
```
