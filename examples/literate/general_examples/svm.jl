# ### Support Vector Machine (SVM)
# We are given two sets of points in ${\bf R}^n$, $\{x_1, \ldots, x_N\}$ and $\{y_1, \ldots, y_M\}$, and wish to find a function $f(x) = w^T x - b$ that linearly separates the points, i.e. $f(x_i) \geq 1$ for $i = 1, \ldots, N$ and $f(y_i) \leq -1$ for $i = 1, \ldots, M$. That is, the points are separated by two hyperplanes, $w^T x - b = 1$ and $w^T x - b = -1$.
#
# Perfect linear separation is not always possible, so we seek to minimize the amount that these inequalities are violated. The violation of point $x_i$ is $\text{max} \{1 + b - w^T x_i, 0\}$, and the violation of point $y_i$ is $\text{max} \{1 - b + w^T y_i, 0\}$. We tradeoff the error $\sum_{i=1}^N \text{max} \{1 + b - w^T x_i, 0\} + \sum_{i=1}^M \text{max} \{1 - b + w^T y_i, 0\}$ with the distance between the two hyperplanes, which we want to be large, via minimizing $\|w\|^2$.
#
# We can write this problem as
# \begin{array}{ll}
#     \mbox{minimize}   & \|w\|^2 + C * (\sum_{i=1}^N \text{max} \{1 + b - w^T x_i, 0\} + \sum_{i=1}^M \text{max} \{1 - b + w^T y_i, 0\}) \\
# \end{array},
# where $w \in {\bf R}^n$ and $b \in {\bf R}$ are our optimization variables.
#
# We can solve the problem as follows.

using Convex
using SCS

#-

## Generate data.
n = 2; # dimensionality of data
C = 10; # inverse regularization parameter in the objective
N = 10; # number of positive examples
M = 10; # number of negative examples

using Distributions
## positive data points
pos = rand(MvNormal([1.0, 2.0], 1.0), N);
## negative data points
neg = rand(MvNormal([-1.0, 2.0], 1.0), M);

#-

function svm(pos, neg, solver=SCSSolver(verbose=0))
    # Create variables for the separating hyperplane w'*x = b.
    w = Variable(n)
    b = Variable()
    # Form the objective.
    obj = sumsquares(w) + C*sum(max(1+b-w'*pos, 0)) + C*sum(max(1-b+w'*neg, 0))
    # Form and solve problem.
    problem = minimize(obj)
    solve!(problem, solver)
    return evaluate(w), evaluate(b)
end;

#-

w, b = svm(pos, neg);

#-

## Plot our results.
using Gadfly
## Generate the separating hyperplane
line_x = -2:0.1:2;
line_y = (-w[1] * line_x + b)/w[2];
## Plot the positive points, negative points, and separating hyperplane.
plot(Scale.y_continuous(minvalue=-1, maxvalue=4),
    layer(x=pos[1,:], y=pos[2,:], Geom.point,
        Theme(default_color=color("green"))),
    layer(x=neg[1,:], y=neg[2,:], Geom.point,
        Theme(default_color=color("red"))),
    layer(x=line_x, y=line_y, Geom.line,
        Theme(default_color=color("blue"))))

