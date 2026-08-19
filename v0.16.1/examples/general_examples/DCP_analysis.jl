# # DCP analysis
using Convex
x = Variable();
y = Variable();
expr = quadoverlin(x - y, 1 - max(x, y))

# We can see from the printing of the expression that this `quadoverlin` (`qol`) atom
# is convex with positive sign. We can query these programmatically using the `vexity`
# and `sign` functions:
println("expression convexity = ", vexity(expr));
println("expression sign = ", sign(expr));
