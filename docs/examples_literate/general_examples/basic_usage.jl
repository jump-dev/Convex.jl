# # Basic Usage

using Convex
using LinearAlgebra
using SCS

# ### Linear program
#
# $$
# \begin{array}{ll}
#   \text{maximize} & c^T x \\
#   \text{subject to} & A x \leq b\\
#   & x \geq 1 \\
#   & x \leq 10 \\
#   & x_2 \leq 5 \\
#   & x_1 + x_4 - x_2 \leq 10 \\
# \end{array}
# $$
#

x = Variable(4)
c = [1; 2; 3; 4]
A = I(4)
b = [10; 10; 10; 10]
p = minimize(dot(c, x)) # or c' * x
p.constraints += A * x <= b
p.constraints += [x >= 1; x <= 10; x[2] <= 5; x[1] + x[4] - x[2] <= 10]
solve!(p, SCS.Optimizer; silent_solver = true)

println(round(p.optval, digits = 2))
println(round.(evaluate(x), digits = 2))
println(evaluate(x[1] + x[4] - x[2]))

# ### Matrix Variables and promotions
#
# $$
# \begin{array}{ll}
#   \text{minimize} & \| X \|_F + y \\
#   \text{subject to} & 2 X \leq 1\\
#   & X' + y \geq 1 \\
#   & X \geq 0 \\
#   & y \geq 0 \\
# \end{array}
# $$
#

X = Variable(2, 2)
y = Variable()
## X is a 2 x 2 variable, and y is scalar. X' + y promotes y to a 2 x 2 variable before adding them
p = minimize(norm(X) + y, 2 * X <= 1, X' + y >= 1, X >= 0, y >= 0)
solve!(p, SCS.Optimizer; silent_solver = true)
println(round.(evaluate(X), digits = 2))
println(evaluate(y))
p.optval

# ### Norm, exponential and geometric mean
#
# $$
# \begin{array}{ll}
#   \text{satisfy} & \| x \|_2 \leq 100 \\
#   & e^{x_1} \leq 5 \\
#   & x_2 \geq 7 \\
#   & \sqrt{x_3 x_4} \geq x_2
# \end{array}
# $$
#

x = Variable(4)
p = satisfy(
    norm(x) <= 100,
    exp(x[1]) <= 5,
    x[2] >= 7,
    geomean(x[3], x[4]) >= x[2],
)
solve!(p, SCS.Optimizer; silent_solver = true)
println(p.status)
evaluate(x)

# ### SDP cone and Eigenvalues

y = Semidefinite(2)
p = maximize(eigmin(y), tr(y) <= 6)
solve!(p, SCS.Optimizer; silent_solver = true)
p.optval

#-

x = Variable()
y = Variable((2, 2))
## SDP constraints
p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
solve!(p, SCS.Optimizer; silent_solver = true)
evaluate(y)

# ### Mixed integer program
#
# $$
# \begin{array}{ll}
#   \text{minimize} & \sum_{i=1}^n x_i \\
#     \text{subject to} & x \in \mathbb{Z}^n \\
#   & x \geq 0.5 \\
# \end{array}
# $$
#

using GLPK
x = Variable(4, :Int)
p = minimize(sum(x), x >= 0.5)
solve!(p, GLPK.Optimizer; silent_solver = true)
evaluate(x)
#-
