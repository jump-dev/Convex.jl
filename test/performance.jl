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
