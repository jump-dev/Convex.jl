@everywhere require("src/CVX.jl")
@everywhere using CVX

TOLERANCE = .0001

# Test 1
x = Variable(1)
p = Problem(:minimize, -x, [x <= 0])
solve!(p)
@assert abs(p.optval - 0) <= TOLERANCE

# Test 2
x = Variable(1)
p = Problem(:minimize, x, [x >= 2, x <= 4])
solve!(p)
@assert abs(p.optval - 2) <= TOLERANCE

# Test 3
x = Variable(1)
p = Problem(:minimize, 2.0 * x, [x >= 2, x <= 4])
solve!(p)
@assert abs(p.optval - 4) <= TOLERANCE

# Test 4
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

# Test 6
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

# Test 9
y = Variable(1)
x = Variable(3)
z = [1.0, 2.0, 3.0] * y
k = -y * [1.0, 2.0, 3.0]
c = [y <= 3.0, y >= 0.0, x >= ones(3), k <= x, x <= z]
o = 3 * y
p = Problem(:minimize,o,c)
solve!(p)
@assert abs(p.optval - 3) < TOLERANCE

# Test 10
X = Variable(2, 2)
c = ones(2, 1)
p = Problem(:minimize, c' * X * c, [X >= ones(2, 2)])
solve!(p)
@assert abs(p.optval - 4) < TOLERANCE

# Test 11
p = Problem(:maximize, c' * X * c, [X <= [1 2; 3 4]])
solve!(p)
@assert abs(p.optval - 10) < TOLERANCE

# Test 12
X = Variable(2, 2)
I = Constant(eye(2))
c = ones(2, 1)
p = Problem(:maximize, c' * (eye(2) + X + I + 1) * c, [X + ones(2, 2) <= [1 2; 3 4]])
solve!(p)
@assert abs(p.optval - 14) < TOLERANCE

# Test 13
x = Variable(1)
p = Problem(:minimize, x, [x >= eye(2)])
solve!(p)
@assert abs(p.optval - 1) < TOLERANCE

# Test 14
x = Variable(1)
p = Problem(:minimize, x, [x + eye(2) >= eye(2)])
solve!(p)
@assert abs(p.optval - 0) < TOLERANCE

# Test 15
x = Variable(1)
c = ones(2, 1)
p = Problem(:minimize, c' * (x + eye(2)) * c, [x + eye(3) >= 2*eye(3), -eye(4) < x])
solve!(p)
@assert abs(p.optval - 6) < TOLERANCE

# Test 16
N = 20
x = Variable(1)
y = Variable(N, N)
c = ones(N, 1)
p = Problem(:minimize, c' * (y + x) * c, [x >= 3, 2y >= 0, y <= x])
solve!(p)
@assert abs(p.optval - 1200) < TOLERANCE

# Test 17
x = Variable(2)
c = ones(2, 1)
p = Problem(:minimize, x' * c, x >= 1)
solve!(p)
@assert abs(p.optval - 2) < TOLERANCE

# Test 18
rows = 2
cols = 3
r = rand(rows, cols)
r_2 = rand(cols, rows)
x = Variable(rows, cols)
c = ones(1, cols)
d = ones(rows, 1)
p = Problem(:minimize, c * x' * d + d' * x * c' + (c * x''''' * d)',
            [x' >= r_2, x >= r, x''' >= r_2, x'' >= r]);
solve!(p)
s = Base.sum(Base.max(r, r_2')) * 3
@assert abs(p.optval - s) < TOLERANCE

# Test 19
rows = 6
cols = 8
n = 2
X = Variable(rows, cols)
A = randn(rows, cols)
c = rand(1, n)
p = Problem(:minimize, c * X[1:n, 5:5+n-1]' * c', X >= A)
solve!(p);
s = c * A[1:n, 5:5+n-1]' * c'
@assert abs(p.optval - s[1]) < TOLERANCE

# Test 20
x = Variable(1, 10)
p = Problem(:minimize, sum(x[2:5]), x >= [1 2 3 4 5 6 7 8 9 10])
solve!(p)
@assert abs(p.optval - 14) < TOLERANCE

# Test 21
x = Variable(10)
a = rand(10, 1)
p = Problem(:maximize, sum(x[2:6]), x <= a)
solve!(p)
@assert abs(p.optval - sum(a[2:6])) < TOLERANCE

# Test 22
x = Variable(10)
a = rand(10, 1)
p = Problem(:minimize, max(x), x >= a)
solve!(p)
@assert abs(p.optval - Base.maximum(a)) < TOLERANCE

# Test 23
x = Variable(10, 10)
y = Variable(10, 10)
a = rand(10, 10)
b = rand(10, 10)
p = Problem(:minimize, max(max(x, y)), [x >= a, y >= b])
solve!(p)
max_a = Base.maximum(a)
max_b = Base.maximum(b)
@assert abs(p.optval - Base.max(max_a, max_b)) < TOLERANCE

# Test 24
x = Variable(1)
a = rand(10, 10)
p = Problem(:maximize, min(x), x <= a)
solve!(p)
@assert abs(p.optval - Base.minimum(a)) < TOLERANCE

# Test 25
x = Variable(10, 10)
y = Variable(10, 10)
a = rand(10, 10)
b = rand(10, 10)
p = Problem(:maximize, min(min(x, y)), [x <= a, y <= b])
solve!(p)
min_a = Base.minimum(a)
min_b = Base.minimum(b)
@assert abs(p.optval - Base.min(min_a, min_b)) < TOLERANCE

# Test 26
x = Variable(3)
a = [-2; 1; 2]
p = Problem(:minimize, sum(pos(x)), [x >= a, x <= 2])
solve!(p)
@assert abs(p.optval - 3) < TOLERANCE

# Test 27
x = Variable(3)
p = Problem(:minimize, norm_inf(x), [-2 <= x, x <= 1])
solve!(p)
@assert abs(p.optval) < TOLERANCE

# Test 28
x1 = Variable(1)
x2 = Variable(1)
p = Problem(:minimize, 4x1 + x2,
            [3x1 + x2 == 3, 4x1 + 3x2 >= 6, x1 + 2x2 <=3, x1 >=0, x2 >=0])
solve!(p)
@assert abs(p.optval - 3.6) < TOLERANCE

# Test 29
x = Variable(4, 4)
y = Variable(4, 6)
p = Problem(:maximize, sum(x) + sum(y), [hcat(x, y) <= 2])
solve!(p)
@assert (p.optval - 80) < TOLERANCE

# Test 30
x = Variable(4, 4)
y = Variable(4, 6)
p = maximize(sum(x) + sum(y), [vertcat(x, y') <= 2])
solve!(p)
@assert (p.optval - 80) < TOLERANCE

# Test
# x = Variable(3)
# p = Problem(:minimize, sum(abs(x)), [-2 <= x, x <= 1])
# solve!(p)
# @assert abs(p.optval) < TOLERANCE

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
