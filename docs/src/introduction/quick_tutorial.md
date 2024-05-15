# Quick Tutorial

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

This problem can be solved in Convex.jl as follows:
```@repl
using Convex, SCS
m = 4;  n = 5
A = randn(m, n); b = randn(m)
x = Variable(n)
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, SCS.Optimizer; silent_solver = true)
problem.status
problem.optval
x.value
```
