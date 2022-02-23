# #  Portfolio Optimization
#
# In this problem, we will find the portfolio allocation that minimizes risk while achieving a given expected return $R_\text{target}$.
#
# Suppose that we know the mean returns $\mu \in \mathbf{R}^n$ and the covariance $\Sigma \in \mathbf{R}^{n \times n}$ of the $n$ assets. We would like to find a portfolio allocation $w \in \mathbf{R}^n$, $\sum_i w_i = 1$, minimizing the *risk* of the portfolio, which we measure as the variance $w^T \Sigma w$ of the portfolio. The requirement that the portfolio allocation achieve the target expected return can be expressed as $w^T \mu >= R_\text{target}$. We suppose further that our portfolio allocation must comply with some lower and upper bounds on the allocation, $w_\text{lower} \leq w \leq w_\text{upper}$.
#
# This problem can be written as
#
# $$
# \begin{array}{ll}
#     \text{minimize}   & w^T \Sigma w \\
#     \text{subject to} & w^T \mu >= R_\text{target} \\
#                       & \sum_i w_i = 1 \\
#                       & w_\text{lower} \leq w \leq w_\text{upper}
# \end{array}
# $$
#
# where $w \in \mathbf{R}^n$ is our optimization variable.

using Convex, SCS

## generate problem data
μ = [11.5; 9.5; 6] / 100          #expected returns
Σ = [
    166 34 58              #covariance matrix
    34 64 4
    58 4 100
] / 100^2

n = length(μ)                   #number of assets

R_target = 0.1
w_lower = 0
w_upper = 0.5;

# If you want to try the optimization with more assets, uncomment and run the next cell. It creates a vector or average returns and a variance-covariance matrix that have scales similar to the numbers above.

#=
using Random
Random.seed!(123)

n = 15                                      #number of assets, CHANGE IT?

μ = (6 .+ (11.5-6)*rand(n))/100             #mean
A = randn(n,n)
Σ = (A * A' + diagm(0=>rand(n)))/500;       #covariance matrix
=#

#-

w = Variable(n)
ret = dot(w, μ)
risk = quadform(w, Σ)

p = minimize(risk, ret >= R_target, sum(w) == 1, w_lower <= w, w <= w_upper)

solve!(p, SCS.Optimizer)

#-

# Optimal portfolio weights:
evaluate(w)

#-

sum(evaluate(w))
