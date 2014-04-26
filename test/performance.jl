function perfTest(N=100000)
  start = time()
  x = Variable(N);
  c = ones(N, 1);
  p = Problem(:minimize, dot(c,x), [c<=x, x<=2*c]);
  println(time() - start)
  solve!(p);
  println(time() - start)
end
