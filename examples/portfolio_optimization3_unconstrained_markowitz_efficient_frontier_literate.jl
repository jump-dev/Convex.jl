# #  Portfolio Optimization - Unconstrained Markowitz Efficient Frontier 
#
# In this problem, we will find the unconstrained portfolio allocation where we introduce the weighting parameter $\lambda(0 \leq \lambda \leq$ 1)and minimize the $\lambda * risk - (1-\lambda)* return$. By varying the values of $\lambda$, we trace out the efficient frontier.  
#
# Suppose that we know the mean returns $R \in \mathbf{R}^n$ of each asset and the covariance $Q \in \mathbf{R}^{n \times n}$ between the assets. Our objective is to find a portfolio allocation that minimizes the *risk* (which we measure as the variance $w^T Q w$) and maximizes the *return* ($w^T R$) of the portfolio of the simulataneously. We suppose further that our portfolio allocation must comply with some lower and upper bounds on the allocation, $w_\mbox{lower} \leq w \leq w_\mbox{upper}$ and also $w \in \mathbf{R}^n$ $\sum_i w_i = 1$.
#
# This problem can be written as
#
# \begin{array}{ll}
#     \mbox{minimize}   & \lambda*w^T Q w - (1-\lambda)*w^T R \\
#     \mbox{subject to} & \sum_i w_i = 1 \\
#                       & w_\mbox{lower} \leq w \leq w_\mbox{upper}
# \end{array}
#
# where $w \in \mathbf{R}^n$ is the vector containing weights allocated to each asset in the efficient frontier.
#
# We can solve this problem as follows.

using Convex, ECOS    #We are using ECOS solver. Install using Pkg.add("ECOS")

## generate problem data
srand(0);     #Set the seed
n = 5;        # Assume that we have portfolio of 5 assets.
R = 5 * randn(n);
A = randn(n, 5);
Q = A * A' + diagm(rand(n));
w_lower = 0;
w_upper = 1;


risk = zeros(2000);   # Initialized the risk and the return vectors.   
ret = zeros(2000);   # lambda varies in the interval(0,1) in the steps of 1/2000.

w = Variable(length(R));

#Defining constraints
c1 = sum(w) == 1;
c2 = w_lower <= w; 
c3 = w <= w_upper;
for i in 1:2000
    λ = i/2000;

    #Defining Objective function
    objective = λ * quadform(w,Q) - (1-λ)* w' *R;
    p = minimize(objective, c1,c2,c3);
    solve!(p, ECOSSolver(verbose = false));
    risk[i] = (w.value' * Q * w.value)[1];
    ret[i] = (w.value'R)[1];
    #println("$i ","$(λ*risk[i] - (1-λ)*ret[i]) ","$p.optval");
    end

using PyPlot            #Install PyPlot if you don't have it installed. Pkg.add("PyPlot")
plot(risk,ret)
title("Markowitz Efficient Frontier");
xlabel("Expected Risk-Variance");
ylabel("Expected Return");

# <img src="efficient_frontier.png">

