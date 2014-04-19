@everywhere require("src/CVX.jl")
@everywhere using CVX


# Test 1
x = Variable(1)
p = Problem(:minimize, -x, [x <= 0])
solve!(p)

# Test 2
x = Variable(1)
p = Problem(:minimize, x, [x >= 2, x <= 4])
solve!(p)

# Test 3
x = Variable(1)
p = Problem(:minimize, 2.0 * x, [x >= 2, x <= 4])
solve!(p)

# Test 4
x = Variable(2)
p = Problem(:minimize, dot([2.0; 2.0], x), [x >= [1.1; 1.1]])
solve!(p)

# Test 5
x = Variable(2)
A = 1.5 * eye(2)
p = Problem(:minimize, dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
solve!(p)

# Test 6
x = Variable(1)
y = Variable(1)
p = Problem(:minimize, x + y, [x >= 3, y >= 2])
solve!(p)

# Test 7
x = Variable(1)
y = Variable(1)
p = Problem(:minimize, 2.0*x - 5.0*y, [100.0 <= x, x >= 200.0, 80.0 <= y, y <= 170.0, y >= -x])
solve!(p)

# Test 8
x = Variable(2)
y = Variable(2)
p = Problem(:minimize, dot([2.0; 2.0], x) + dot([2.0; 2.0], y), [x >= [1.1; 1.1], y >= [1.1; 1.1]])
solve!(p)

# TODO: take care of this
# p = Problem(:minimize, 2*x - 5*y, [100 <= x, x >= 200, 80 <= y, y <= 170, y >= -x + 200])
# solve!(p)

# sol = ecos_solve(n=2, m=2, p=0, G=[-1 -2; -2 -1], c=[1.0, 1.0], h=[-2.0, -2.0])
#TOL = .0001

#x = Parameter(3)
#y = Variable(3,4)
#z = x + 2
#w = 3*y
#q = rand(3,3)*w + y

#y = Variable()
#x = Variable(3)
#z = y*[1.0,2.0,3.0]
#c = [y <= 3.0, y >= 0.0, x >= ones(3), x <= z]
#o = 3*y
#p = Problem(:minimize,o,c)
#solve!(p)
#@assert abs(p.optval - 3) < TOL
#println(x.value)
#@assert x.value - y.value*[1.0,2.0,3.0] <= 0

# try some atoms
#y = Variable()
#p = Problem(:minimize,abs(y),[y >= 1])
#@assert abs(solve!(p) - 1) <= TOL

# x = Variable()
# p = Problem(:minimize,kl_div(x,y),[y >= 1 , x <= 0])
# print
# @assert abs(solve!(p) - 1) <= TOL
