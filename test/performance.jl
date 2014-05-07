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

# Assuming speed of calling rand, sprand, etc is negligible
function test2(m=1e4, n=1e4, q=1e4)
  x = Variable(n);
  b = Variable(m);
  d = Variable(q);
  A = rand(m,n);
  C = rand(q,n);
  c = rand(n);
  p = Problem(:minimize, c' * x, [A*x <=b, C*x <=d])
  solve!(p)
  println(p.optval)
  return p
end

function test2sparse(m=1e5, n=1e5, q=1e5, p=0.1)
  x = Variable(n);
  b = Variable(m);
  d = Variable(q);
  A = sprand(m,n,p);
  C = sprand(q,n,p);
  c = rand(n);
  p = Problem(:minimize, c' * x, [A*x <=b, C*x <=d])
  solve!(p)
  println(p.optval)
  return p
end

# http://see.stanford.edu/materials/lsocoee364a/hw4sol.pdf
function test3(m=1e4,n=1e4) 
  t = Variable();
  A = rand(m,n);
  x = rand(n); # TODO: The more interesting problem is to make x also a variable
  # but then we don't have the solution analytically
  # we can compare it to minimize norm(A*x-b, 'inf') in Matlab with variable x
  b = rand(m);
  p = Problem(:minimize, t, [A*x - b <= t*ones(m,1), -t*ones(m,1) <= A*x - b])
  solve!(p)
  println(p.optval)
  soln = norm(A*x - b, Inf)
  println(soln)
  println(p.optval - soln)
  return p
end

# http://see.stanford.edu/materials/lsocoee364a/hw4sol.pdf
function test4(m=1e4, n=1e4)
  s = Variable(m);
  x = rand(n); # TODO: This is also more interesting with x = Variable(n)
  # We can compare it to minimize norm(A*x-b, 1) in Matlab with variable x
  A = rand(m,n);
  b = rand(m);
  p = Problem(:minimize, sum(s), [A*x - b <= s, A*x - b >= -s])
  solve!(p)
  return p
end