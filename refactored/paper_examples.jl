using ECOS, SCS

# Summation.
x = Variable();
e = 0;
@time begin
  for i = 1:1000
    e += x;
  end
  p = minimize(e);
end
@time solve!(p, ECOS.ECOSMathProgModel())

# Indexing.
x = Variable(1000, 1);
e = 0;
@time begin
  for i = 1:1000
    e += x[i];
  end
  p = minimize(e);
end
@time solve!(p, ECOS.ECOSMathProgModel())

# Matrix constraints.
n, m, p = 100, 100, 100
X = Variable(m, n)
A = randn(p, m)
b = randn(p, n)
@time begin
  p = minimize(vecnorm(X), A * X == b)
end
@time solve!(p, ECOS.ECOSMathProgModel())

# Transpose.
X = Variable(1000, 1000)
A = randn(1000, 1000)
@time begin
  p = minimize(norm2(X - A), X' == X)
end
@time solve!(p)
