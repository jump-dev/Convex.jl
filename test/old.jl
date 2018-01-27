function spdiags(B,d,m,n)
  d = d[:]
  p = length(d)

  len = zeros(p+1,1)
  for k = 1:p
      len[k+1] = int(len[k]+length(Base.max(1,1-d[k]): Base.min(m,n-d[k])))
  end
  a = zeros(int(len[p+1]),3)
  for k = 1:p
      # Append new d[k]-th diagonal to compact form
      i = Base.max(1,1-d[k]):Base.min(m,n-d[k])
      a[(int(len[k])+1):int(len[k+1]),:] = [i i+d[k] B[i+(m>=n)*d[k],k]]
  end

  A = sparse(int(a[:,1]),int(a[:,2]),a[:,3],m,n)

  return A
end

y = readcsv("snp_500.txt")
n = length(y)

e = ones(n, 1)
D = spdiags([e -2*e e], 0:2, n-2, n)

lambda = 50
start = time()
D*x
println(time() - start)

x = Variable(n)
p = minimize(lambda*sum(D*x))
println(time() - start)
solve!(p)
p.optval
println(time() - start)