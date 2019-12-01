# # Binary (or 0-1) knapsack problem
# Given a knapsack of some capacity $C$ and $n$ objects with object $i$ having weight $w_i$ and profit $p_i$, the goal is to choose some subset of the objects that can fit in the knapsack (i.e. the sum of their weights is no more than $C$) while maximizing profit.
#
# This can be formulated as a mixed-integer program as:
#
# $$
# \begin{array}{ll}
#   \text{maximize} & x' p \\
#     \text{subject to} & x \in \{0, 1\} \\
#   & w' x \leq C \\
# \end{array}
# $$
#
# where $x$ is a vector is size $n$ where $x_i$ is one if we chose to keep the object in the knapsack, 0 otherwise.

## Data taken from http://people.sc.fsu.edu/~jburkardt/datasets/knapsack_01/knapsack_01.html
w = [23; 31; 29; 44; 53; 38; 63; 85; 89; 82]
C = 165 
p =  [92; 57; 49; 68; 60; 43; 67; 84; 87; 72];
n = length(w)

#-

using Convex, GLPK
x = Variable(n, :Bin)
problem = maximize(dot(p, x), dot(w, x) <= C)
solve!(problem, GLPK.Optimizer())
evaluate(x)
