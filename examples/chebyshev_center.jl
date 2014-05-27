# Boyd & Vandenberghe, "Convex Optimization"
# Joëlle Skaf - 08/16/05
#
# Adapted for CVX.jl by Karanveer Mohan and David Zeng - 26/05/14
#
# The goal is to find the largest Euclidean ball (i.e. its center and
# radius) that lies in a polyhedron described by linear inequalites in this
# fashion: P = {x : a_i'*x <= b_i, i=1,...,m} where x is in R^2

# Generate the input data
a1 = [ 2;  1];
a2 = [ 2; -1];
a3 = [-1;  2];
a4 = [-1; -2];
b = ones(4, 1);

# Create and solve the model
r = Variable(1)
x_c = Variable(2)
p = maximize(r)
p.constraints += a1' * x_c + r * Base.norm(a1, 2) <= b[1];
p.constraints += a2' * x_c + r * Base.norm(a2, 2) <= b[2];
p.constraints += a3' * x_c + r * Base.norm(a3, 2) <= b[3];
p.constraints += a4' * x_c + r * Base.norm(a4, 2) <= b[4];
solve!(p)
p.optval

# Generate the figure
x = linspace(-1.5, 1.5);
theta = 0:pi/100:2*pi;
using Gaston
Gaston.set_terminal("x11")
plot(x, -x * a1[1] / a1[2] + b[1] / a1[2], "color", "black",
     x, -x * a2[1]/ a2[2] + b[2] / a2[2], "color", "red",
     x, -x * a3[1]/ a3[2] + b[3] / a3[2], "color", "green",
     x, -x * a4[1]/ a4[2] + b[4] / a4[2], "color", "blue",
     x_c.value[1] + r.value * cos(theta), x_c.value[2] + r.value * sin(theta), "color", "red",
     "title", "Largest Euclidean ball lying in a 2D polyhedron",
    );
