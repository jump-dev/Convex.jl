require("ast.jl")
TOL = .0001

x = Parameter(3)
y = Variable(3,4)
z = x + 2
w = 3*y
q = rand(3,3)*w + y

y = Variable()
x = Variable(3)
z = y*[1.0,2.0,3.0]
c = [y <= 3.0, y >= 0.0, x >= ones(3), x <= z]
o = 3*y
p = Problem(:minimize,o,c)
solve!(p)
@assert abs(p.optval - 3) < TOL
# println(x.value)
#@assert x.value - y.value*[1.0,2.0,3.0] <= 0

# try some atoms
y = Variable()
p = Problem(:minimize,abs(y),[y >= 1])
@assert abs(solve!(p) - 1) <= TOL

x = Variable()
p = Problem(:minimize,kl_div(x,y),[y >= 1 , x <= 0])
@assert abs(solve!(p) - 1) <= TOL