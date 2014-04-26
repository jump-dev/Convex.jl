@everywhere require("src/CVX.jl")
@everywhere using CVX

TOLERANCE = .0001

# Test 1
x = Variable(1)
p = Problem(:minimize, -x, [x <= 0])
solve!(p)
@assert abs(p.optval - 0) <= TOLERANCE

# # Test 2
x = Variable(1)
p = Problem(:minimize, x, [x >= 2, x <= 4])
solve!(p)
@assert abs(p.optval - 2) <= TOLERANCE

# Test 3
x = Variable(1)
p = Problem(:minimize, 2.0 * x, [x >= 2, x <= 4])
solve!(p)
@assert abs(p.optval - 4) <= TOLERANCE

# # Test 4
x = Variable(2)
p = Problem(:minimize, dot([2.0; 2.0], x), [x >= [1.1; 1.1]])
solve!(p)
@assert abs(p.optval - 4.4) <= TOLERANCE


# Test 5
x = Variable(2)
A = 1.5 * eye(2)
p = Problem(:minimize, dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
solve!(p)
@assert abs(p.optval - 2.9333) <= TOLERANCE

# # Test 6
x = Variable(1)
y = Variable(1)
p = Problem(:minimize, x + y, [x >= 3, y >= 2])
solve!(p)
@assert abs(p.optval - 5) <= TOLERANCE

# Test 7
x = Variable(1)
y = Variable(1)
p = Problem(:minimize, 2.0*x - 5.0*y, [100.0 <= x, x <= 200.0, 80.0 <= y, y <= 170.0, y >= -x])
solve!(p)
@assert abs(p.optval + 650) <= TOLERANCE

# Test 8
x = Variable(2)
p = Problem(:minimize, x[1] + x[2], [x >= 1])
solve!(p)
@assert abs(p.optval - 2) <= TOLERANCE

y = Variable(1)
x = Variable(2)
z = y*[1; 2]
c = [y >= 0.0, x >= ones(2), x <= z]
p = Problem(:minimize,y,c)
# solve!(p)

# y = Variable(1)
# x = Variable(3)
# z = y*[1.0,2.0,3.0]
# c = [y <= 3.0, y >= 0.0, x >= ones(3), x <= z]
# o = 3*y
# p = Problem(:minimize,o,c)
# solve!(p)
# @assert abs(p.optval - 3) < TOL
# println(x.value)
# @assert x.value - y.value*[1.0,2.0,3.0] <= 0

# try some atoms
#y = Variable()
#p = Problem(:minimize,abs(y),[y >= 1])
#@assert abs(solve!(p) - 1) <= TOL

# x = Variable()
# p = Problem(:minimize,kl_div(x,y),[y >= 1 , x <= 0])
# print
# @assert abs(solve!(p) - 1) <= TOL
