# #  Portfolio Optimization
#
# In this problem, we will find the portfolio allocation that minimizes risk while achieving a given expected return $R_\mbox{target}$.
#
# Suppose that we know the mean returns $R \in \mathbf{R}^n$ and the covariance $Q \in \mathbf{R}^{n \times n}$ of the $n$ assets. We would like to find a portfolio allocation $x \in \mathbf{R}^n$, $\sum_i x_i = 1$, minimizing the *risk* of the portfolio, which we measure as the variance $x^T Q x$ of the portfolio. The requirement that the portfolio allocation achieve the target expected return can be expressed as $x^T R >= R_\mbox{target}$. We suppose further that our portfolio allocation must comply with some lower and upper bounds on the allocation, $x_\mbox{lower} \leq x \leq x_\mbox{upper}$.
#
# This problem can be written as
#
# \begin{array}{ll}
#     \mbox{minimize}   & x^T Q x \\
#     \mbox{subject to} & x^T R >= R_\mbox{target} \\
#                       & \sum_i x_i = 1 \\
#                       & x_\mbox{lower} \leq x \leq x_\mbox{upper}
# \end{array}
#
# where $x \in \mathbf{R}^n$ is our optimization variable.
#
# We can solve this problem as follows.

using Convex, ECOS

## generate problem data
srand(0)
n = 50
R = 5*randn(n)
A = randn(n, 5)
Q = A * A' + diagm(rand(n))
R_target = 5
x_lower = 0
x_upper = 1

x = Variable(length(R))
p = minimize(quadform(x, Q), 
             x' * R >= R_target, 
             sum(x) == 1, 
             x_lower <= x, 
             x <= x_upper )

solve!(p, ECOSSolver(verbose = false)) 

## the minimal risk
p.optval

# We see that we can achieve an extremely low risk portfolio (with variance .025) with the desired expected return. 

#-

# The optimal portfolio invests in only about half of the assets.

sum(x.value.>1e-4)

# Let's take a look at the optimal portfolio we chose:

plot(x=1:n,y=x.value,Geom.bar,Guide.xlabel("Asset Index"),Guide.ylabel("Fraction of Portfolio"))

