@everywhere require("src/CVX.jl")
@everywhere using CVX

TOLERANCE = .0001

# # Test 1
# x = Variable(1)
# p = Problem(:minimize, -x, [x <= 0])
# solve!(p)
# @assert abs(p.optval - 0) <= TOLERANCE

# # # Test 2
# x = Variable(1)
# p = Problem(:minimize, x, [x >= 2, x <= 4])
# solve!(p)
# @assert abs(p.optval - 2) <= TOLERANCE

# # Test 3
# x = Variable(1)
# p = Problem(:minimize, 2.0 * x, [x >= 2, x <= 4])
# solve!(p)
# @assert abs(p.optval - 4) <= TOLERANCE

# # Test 4
# x = Variable(2)
# p = Problem(:minimize, dot([2.0; 2.0], x), [x >= [1.1; 1.1]])
# solve!(p)
# @assert abs(p.optval - 4.4) <= TOLERANCE


# # Test 5
# x = Variable(2)
# A = 1.5 * eye(2)
# p = Problem(:minimize, dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
# solve!(p)
# @assert abs(p.optval - 2.9333) <= TOLERANCE

# # Test 6
# x = Variable(1)
# y = Variable(1)
# p = Problem(:minimize, x + y, [x >= 3, y >= 2])
# solve!(p)
# @assert abs(p.optval - 5) <= TOLERANCE

# # Test 7
# x = Variable(1)
# y = Variable(1)
# p = Problem(:minimize, 2.0*x - 5.0*y, [100.0 <= x, x <= 200.0, 80.0 <= y, y <= 170.0, y >= -x])
# solve!(p)
# @assert abs(p.optval + 650) <= TOLERANCE

# # Test 8
# x = Variable(2)
# p = Problem(:minimize, x[1] + x[2], [x >= 1])
# solve!(p)
# @assert abs(p.optval - 2) <= TOLERANCE

# # Test 9
# y = Variable(1)
# x = Variable(3)
# z = [1.0, 2.0, 3.0] * y
# k = -y * [1.0, 2.0, 3.0]
# c = [y <= 3.0, y >= 0.0, x >= ones(3), k <= x, x <= z]
# o = 3 * y
# p = Problem(:minimize,o,c)
# solve!(p)
# @assert abs(p.optval - 3) < TOLERANCE

# # Test 10
# X = Variable(2, 2)
# c = ones(2, 1)
# p = Problem(:minimize, c' * X * c, [X >= ones(2, 2)])
# solve!(p)
# @assert abs(p.optval - 4) < TOLERANCE

# # Test 11
# p = Problem(:maximize, c' * X * c, [X <= [1 2; 3 4]])
# solve!(p)
# @assert abs(p.optval - 10) < TOLERANCE

# # Test 12
# X = Variable(2, 2)
# I = Constant(eye(2))
# c = ones(2, 1)
# p = Problem(:maximize, c' * (eye(2) + X + I + 1) * c, [X + ones(2, 2) <= [1 2; 3 4]])
# solve!(p)
# @assert abs(p.optval - 14) < TOLERANCE

# # Test 13
# x = Variable(1)
# p = Problem(:minimize, x, [x >= eye(2)])
# solve!(p)
# @assert abs(p.optval - 1) < TOLERANCE

# # Test 14
# x = Variable(1)
# p = Problem(:minimize, x, [x + eye(2) >= eye(2)])
# solve!(p)
# @assert abs(p.optval - 0) < TOLERANCE

# Test 15
# x = Variable(1)
# c = ones(2, 1)
# p = Problem(:minimize, c' * (x + eye(2)) * c, [x + eye(3) >= 2*eye(3), -eye(4) < x])
# solve!(p)
# @assert abs(p.optval - 6) < TOLERANCE

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

# x = Variable(1)
# p = Problem(:minimize, x, [eye(2) + x >= ones(2, 2)])
# solve!(p)
# p = Problem(:minimize, y, [X >= 1, y >= X])
# println(x.value)
# @assert x.value - y.value*[1.0,2.0,3.0] <= 0

# try some atoms
#y = Variable()
#p = Problem(:minimize,abs(y),[y >= 1])
#@assert abs(solve!(p) - 1) <= TOLERANCE

# x = Variable()
# p = Problem(:minimize,kl_div(x,y),[y >= 1 , x <= 0])
# print
# @assert abs(solve!(p) - 1) <= TOLERANCE
