# # Optimal advertising

# This example is taken from <https://web.stanford.edu/~boyd/papers/pdf/cvx_applications.pdf>.

# Setup:
#
# * We have $m$ adverts and $n$ time slots
# * The total traffic in time slot $t$ is $T_t$
# * The number of ad $i$ displayed in period $t$ is $D_{it} \geq 0$
# * We require $\sum_{i=1}^m D_{it} \leq T_t$ since we cannot show more than $T_t$ ads during time slot $t$.
# * We require $\sum_{t=1}^n D_{it} \geq c_i$ to fulfill a contract to show advertisement $i$ at least $c_i$ times.
#
# Goal: choose $D_{it}$.
#
# For some empirical $P_{it}$ with $0 \leq P_{it} \leq 1$, we obtain $C_{it} = P_{it}D_{it}$ clicks for ad $i$, which pays us some number $R_i > 0$ up to a budget $B_i$.
# The ad revenue for ad $i$ is $S_i = \min( R_i \sum_t C_{it}, B_i )$ which is concave in $D$.
# We aim to maximize $\sum_i S_i$.

using Random
using Distributions: LogNormal
Random.seed!(1);

m = 5; # number of adverts
n = 24; # number of time slots
SCALE = 10000;
B = rand(LogNormal(8), m) .+ 10000;
B = round.(B, digits = 3); # Budget

P_ad = rand(m);
P_time = rand(1, n);
P = P_ad * P_time;

T = sin.(range(-2 * pi / 2, stop = 2 * pi - 2 * pi / 2, length = n)) * SCALE;
T .+= -minimum(T) + SCALE; # traffic
c = rand(m); # contractual minimum
c *= 0.6 * sum(T) / sum(c);
c = round.(c, digits = 3);
R = [rand(LogNormal(minimum(c) / c[i]), 1) for i in 1:m]; # revenue

#-

## Form and solve the optimal advertising problem.
using Convex, SCS;
D = Variable(m, n);
Si = vcat([min(R[i] * dot(P[i, :], D[i, :]'), B[i]) for i in 1:m]...);
problem =
    maximize(sum(Si), [D >= 0, sum(D, dims = 1)' <= T, sum(D, dims = 2) >= c]);
solve!(problem, SCS.Optimizer; silent = true)

#-

# Plot traffic.
using Plots
plot(1:length(T), T, xlabel = "hour", ylabel = "Traffic")

#-

# Plot P.
heatmap(P)

#-

# Plot optimal D.
heatmap(evaluate(D))
