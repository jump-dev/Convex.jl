# # Paper examples

using Convex, ECOS

# Summation.
println("Summation example")
x = Variable();
e = 0;
@time begin
  for i = 1:1000
    global e
    e += x;
  end
  p = minimize(e, x>=1);
end
@time solve!(p, ECOSSolver())

# Indexing.
println("Indexing example")
x = Variable(1000, 1);
e = 0;
@time begin
  for i = 1:1000
    global e
    e += x[i];
  end
  p = minimize(e, x >= ones(1000, 1));
end
@time solve!(p, ECOSSolver())

# Matrix constraints.
println("Matrix constraint example")
n, m, p = 100, 100, 100
X = Variable(m, n);
A = randn(p, m);
b = randn(p, n);
@time begin
  p = minimize(norm(vec(X)), A * X == b);
end
@time solve!(p, ECOSSolver())

# Transpose.
println("Transpose example")
X = Variable(5, 5);
A = randn(5, 5);
@time begin
  p = minimize(norm2(X - A), X' == X);
end
@time solve!(p, ECOSSolver())

n = 3
A = randn(n, n);
#@time begin
  X = Variable(n, n);
  p = minimize(vecnorm(X' - A), X[1,1] == 1);
  solve!(p, ECOSSolver())
#end
