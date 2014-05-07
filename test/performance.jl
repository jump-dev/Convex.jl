import Base.norm

function perfTest(N=100000)
  start = time()
  x = Variable(N);
  c = ones(N, 1);
  p = Problem(:minimize, dot(c,x), [c<=x, x<=2*c]);
  println(time() - start)
  solve!(p);
  println(time() - start)
end

#TODO code so inefficient even python is faster
function test()
  start = time()
  N = 200
  x = Variable(1)
  y = Variable(N, N)
  c = ones(N, 1)
  p = Problem(:minimize, c' * (y + x) * c, [x >= 3, 2*y>=0, y <= x])
  println(time() - start)
  solve!(p)
  println(p.optval)
  println(time() - start)
  return p
end

# Assuming speed of calling randn, sprandn, etc is negligible for now...
# This is untrue. We'll soon generate data files.
# Also we'll need lucky data for this to be feasible.
# test2(10, 15, 10); sometimes gets feasible, bounded problems
function test2(m=1000, n=1000, q=1000)
  x = Variable(n);
  A = randn(m,n);
  b = randn(m);
  C = randn(q,n);
  d = C * randn(n);
  c = randn(n);
  p = Problem(:minimize, c' * x, [A*x <=b, C*x == d]);
  solve!(p);
  println(p.optval)
  return p
end

# test2sparse(10,10,5,0.4); sometimes gets bounded, feasible problems
function test2sparse(m=1000, n=1000, q=1000, p=0.1)
  x = Variable(n);
  A = sprandn(m,n,p);
  b = randn(m);
  C = sprandn(q,n,p);
  d = C * randn(n);
  c = randn(n);
  p = Problem(:minimize, c' * x, [A*x <=b, C*x ==d]);
  solve!(p);
  println(p.optval)
  return p
end

# http://see.stanford.edu/materials/lsocoee364a/hw4sol.pdf
# close to optimal with 10, 10
function test3(m=1000,n=1000) 
  t = Variable(1);
  A = randn(m,n);
  x = randn(n); # TODO: The more interesting problem is to make x also a variable
  # but then we don't have the solution analytically
  # we can compare it to minimize norm(A*x-b, 'inf') in Matlab with variable x
  b = randn(m);
  p = Problem(:minimize, t, [A*x - b <= t*ones(m,1), -t*ones(m,1) <= A*x - b]);
  solve!(p);
  println(p.optval)
  soln = Base.norm(A*x - b, Inf)
  println("Solution and deviation from solution")
  println(soln)
  println(p.optval - soln)
  return p
end

# http://see.stanford.edu/materials/lsocoee364a/hw4sol.pdf
# we get close to optimal with 10, 10
function test4(m=1000, n=1000)
  s = Variable(m);
  x = randn(n); # TODO: This is also more interesting with x = Variable(n)
  # We can compare it to minimize norm(A*x-b, 1) in Matlab with variable x
  A = randn(m,n);
  b = randn(m);
  p = Problem(:minimize, sum(s), [A*x - b <= s, A*x - b >= -s])
  solve!(p)
  return p
end