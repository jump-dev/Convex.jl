# # Support vector machine
# ### Support Vector Machine (SVM)
# We are given two sets of points in ${\bf R}^n$, $\{x_1, \ldots, x_N\}$ and $\{y_1, \ldots, y_M\}$, and wish to find a function $f(x) = w^T x - b$ that linearly separates the points, i.e. $f(x_i) \geq 1$ for $i = 1, \ldots, N$ and $f(y_i) \leq -1$ for $i = 1, \ldots, M$. That is, the points are separated by two hyperplanes, $w^T x - b = 1$ and $w^T x - b = -1$.
#
# Perfect linear separation is not always possible, so we seek to minimize the amount that these inequalities are violated. The violation of point $x_i$ is $\text{max} \{1 + b - w^T x_i, 0\}$, and the violation of point $y_i$ is $\text{max} \{1 - b + w^T y_i, 0\}$. We tradeoff the error $\sum_{i=1}^N \text{max} \{1 + b - w^T x_i, 0\} + \sum_{i=1}^M \text{max} \{1 - b + w^T y_i, 0\}$ with the distance between the two hyperplanes, which we want to be large, via minimizing $\|w\|^2$.
#
# We can write this problem as
#
# $$
# \begin{array}{ll}
#     \text{minimize}   & \|w\|^2 + C * (\sum_{i=1}^N \text{max} \{1 + b - w^T x_i, 0\} + \sum_{i=1}^M \text{max} \{1 - b + w^T y_i, 0\})
# \end{array},
# $$
#
# where $w \in {\bf R}^n$ and $b \in {\bf R}$ are our optimization variables.
#
# We can solve the problem as follows.

using Convex, SCS

#-

## Generate data.
n = 2; # dimensionality of data
C = 10; # inverse regularization parameter in the objective
N = 10; # number of positive examples
M = 10; # number of negative examples

using Distributions: MvNormal
## positive data points
pos_data = rand(MvNormal([1.0, 2.0], 1.0), N);
## negative data points
neg_data = rand(MvNormal([-1.0, 2.0], 1.0), M);

#-

function svm(pos_data, neg_data, solver=SCS.Optimizer(verbose=0))
    ## Create variables for the separating hyperplane w'*x = b.
    w = Variable(n)
    b = Variable()
    ## Form the objective.
    obj = sumsquares(w) + C*sum(max(1+b-w'*pos_data, 0)) + C*sum(max(1-b+w'*neg_data, 0))
    ## Form and solve problem.
    problem = minimize(obj)
    solve!(problem, solver)
    return evaluate(w), evaluate(b)
end;

#-

w, b = svm(pos_data, neg_data);

#-

## Plot our results.
using Plots
# Generate the separating hyperplane
line_x = -2:0.1:2;
line_y = (-w[1] * line_x .+ b)/w[2];
# Plot the positive points, negative points, and separating hyperplane.
plot(pos_data[1,:], pos_data[2,:], st=:scatter, label="Positive points")
plot!(neg_data[1,:], neg_data[2,:], st=:scatter, label="Negative points")
plot!(line_x, line_y, label="Separating hyperplane")
