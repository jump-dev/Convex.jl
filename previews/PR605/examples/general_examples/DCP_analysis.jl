# # DCP analysis
using Convex
x = Variable();
y = Variable();
expr = quadoverlin(x - y, 1 - max(x, y));
println("expression curvature = ", vexity(expr));
println("expression sign = ", sign(expr));
