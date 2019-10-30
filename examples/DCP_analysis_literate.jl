## DCP analysis in Convex.jl
using Convex
x = Variable();
y = Variable();
expr = quad_over_lin(x - y, 1 - max(x, y));
println("expression curvature = ", vexity(expr));
println("expression sign = ", sign(expr));

