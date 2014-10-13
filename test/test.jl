using Convex

TOLERANCE = .001

# Test 1
x = Variable(1)
p = minimize(-x, [x <= 0])
solve!(p)
@assert abs(p.optval - 0) <= TOLERANCE

# Test 2
x = Variable(1)
p = minimize(x, [x >= 2, x <= 4])
solve!(p)
@assert abs(p.optval - 2) <= TOLERANCE

# Test 3
x = Variable(1)
p = minimize(2.0 * x, [x >= 2, x <= 4])
solve!(p)
@assert abs(p.optval - 4) <= TOLERANCE

# Test 4
x = Variable(2)
constr = x >= [1.1; 1.1]
p = minimize(dot([2.0; 2.0], x), constr)
solve!(p)
@assert abs(p.optval - 4.4) <= TOLERANCE

# Test 5
x = Variable(2)
A = 1.5 * eye(2)
p = minimize(dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
solve!(p)
@assert abs(p.optval - 2.9333) <= TOLERANCE

# Test 6
x = Variable(1)
y = Variable(1)
p = minimize(x + y, [x >= 3, y >= 2])
solve!(p)
@assert abs(p.optval - 5) <= TOLERANCE

# Test 7
x = Variable(1)
y = Variable(1)
p = minimize(2.0*x - 5.0*y, [100.0 <= x, x <= 200.0, 80.0 <= y, y <= 170.0, y >= -x])
solve!(p)
@assert abs(p.optval + 650) <= TOLERANCE

# Test 8
x = Variable(2)
p = minimize(x[1] + x[2], [x >= 1])
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
p = minimize(c' * X * c, [X >= ones(2, 2)])
solve!(p)
@assert abs(p.optval - 4) < TOLERANCE

# Test 11
p = maximize(c' * X * c, [X <= [1 2; 3 4]])
solve!(p)
@assert abs(p.optval - 10) < TOLERANCE

# Test 12
X = Variable(2, 2)
I = Constant(eye(2))
c = ones(2, 1)
p = maximize(c' * (eye(2) + X + I .+ 1) * c, [X + ones(2, 2) <= [1 2; 3 4]])
solve!(p)
@assert abs(p.optval - 14) < TOLERANCE

# Test 13
x = Variable(1)
p = minimize(x, [x >= eye(2)])
solve!(p)
@assert abs(p.optval - 1) < TOLERANCE

# Test 14
x = Variable(1)
p = minimize(x, [x .+ eye(2) >= eye(2)])
solve!(p)
@assert abs(p.optval - 0) < TOLERANCE

# Test 15
x = Variable(1)
c = ones(2, 1)
p = minimize(c' * (x .+ eye(2)) * c, [x .+ eye(3) >= 2 * eye(3), -eye(4) <= x])
solve!(p)
@assert abs(p.optval - 6) < TOLERANCE

# Test 16
N = 20
x = Variable(1)
y = Variable(N, N)
c = ones(N, 1)
p = minimize(c' * (y .+ x) * c, [x >= 3, 2y >= 0, y <= x])
solve!(p)
@assert abs(p.optval - 1200) < TOLERANCE

# Test 17
x = Variable(2)
c = ones(2, 1)
p = minimize(x' * c, x >= 1)
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
p = minimize(c * x' * d + d' * x * c' + (c * x''''' * d)',
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
p = minimize(c * X[1:n, 5:5+n-1]' * c', X >= A)
solve!(p);
s = c * A[1:n, 5:5+n-1]' * c'
@assert abs(p.optval - s[1]) < TOLERANCE

# Test 20
x = Variable(1, 10)
p = minimize(sum(x[2:5]), x >= [1 2 3 4 5 6 7 8 9 10])
solve!(p)
@assert abs(p.optval - 14) < TOLERANCE

# Test 21
x = Variable(10)
a = rand(10, 1)
p = maximize(sum(x[2:6]), x <= a)
solve!(p)
@assert abs(p.optval - sum(a[2:6])) < TOLERANCE

# Test 22
x = Variable(10)
a = rand(10, 1)
p = minimize(maximum(x), x >= a)
solve!(p)
@assert abs(p.optval - Base.maximum(a)) < TOLERANCE

# Test 23
x = Variable(10, 10)
y = Variable(10, 10)
a = rand(10, 10)
b = rand(10, 10)
p = minimize(maximum(max(x, y)), [x >= a, y >= b])
solve!(p)
max_a = Base.maximum(a)
max_b = Base.maximum(b)
@assert abs(p.optval - Base.max(max_a, max_b)) < TOLERANCE

# Test 24
x = Variable(1)
a = rand(10, 10)
p = maximize(minimum(x), x <= a)
solve!(p)
@assert abs(p.optval - Base.minimum(a)) < TOLERANCE

# Test 25
x = Variable(10, 10)
y = Variable(10, 10)
a = rand(10, 10)
b = rand(10, 10)
p = maximize(minimum(min(x, y)), [x <= a, y <= b])
solve!(p)
min_a = Base.minimum(a)
min_b = Base.minimum(b)
@assert abs(p.optval - Base.min(min_a, min_b)) < TOLERANCE

# Test 26
x = Variable(3)
a = [-2; 1; 2]
p = minimize(sum(pos(x)), [x >= a, x <= 2])
solve!(p)
@assert abs(p.optval - 3) < TOLERANCE

# Test 27
x = Variable(3)
p = minimize(norm_inf(x), [-2 <= x, x <= 1])
solve!(p)
@assert abs(p.optval) < TOLERANCE

# Test 28
x1 = Variable(1)
x2 = Variable(1)
p = minimize(4x1 + x2,
            [3x1 + x2 == 3, 4x1 + 3x2 >= 6, x1 + 2x2 <=3, x1 >=0, x2 >=0])
solve!(p)
@assert abs(p.optval - 3.6) < TOLERANCE

# Test 29
x = Variable(4, 4)
y = Variable(4, 6)
p = maximize(sum(x) + sum(y), [hcat(x, y) <= 2])
solve!(p)
@assert abs(p.optval - 80) < TOLERANCE

# Test 30
x = Variable(4, 4)
y = Variable(4, 6)
p = maximize(sum(x) + sum(y), [vertcat(x, y') <= 2])
solve!(p)
@assert abs(p.optval - 80) < TOLERANCE

# Test 31
x = Variable(4, 4)
y = Variable(4, 6)
z = Variable(1)
c = ones(4, 1)
d = 2 * ones(6, 1)
constraints = [hcat(x, y) <= 2, z <= 0, z <= x, 2z >= -1]
objective = sum(x .+ z) + min(y) + c' * y * d
p = maximize(objective, constraints)
solve!(p)
@assert abs(p.optval - 130) < TOLERANCE

# Test 32
x = Variable(2, 1)
A = [1 2; 2 1; 3 4]
b = [2; 3; 4]
p = minimize(norm_2(A * x + b))
solve!(p)
@assert abs(p.optval - 0.64888) < TOLERANCE

# Test 33
x = Variable(2, 1)
A = [1 2; 2 1; 3 4]
b = [2; 3; 4]
lambda = 1
p = minimize(norm_2(A * x + b) + lambda * norm_2(x), x >= 1)
solve!(p)
@assert abs(p.optval - 14.9049) < TOLERANCE

# Test 34
x = Variable(2, 1)
A = [1 2; 2 1; 3 4]
b = [2; 3; 4]
p = minimize(sum_squares(A*x + b))
solve!(p)
@assert abs(p.optval - 0.42105) < TOLERANCE

# Test 35
x = Variable(3)
p = minimize(sum(abs(x)), [-2 <= x, x <= 1]);
solve!(p)
@assert abs(p.optval) < TOLERANCE

# Test 36
x = Variable(2, 1)
A = [1 2; 2 1; 3 4]
b = [2; 3; 4]
lambda = 1
p = minimize(norm_2(A * x + b) + lambda * norm_1(x), x >= 1)
solve!(p)
@assert abs(p.optval - 15.4907) < TOLERANCE

# Test 37
x = Variable(3, 1)
A = [0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
p = minimize(quad_form(x, A), [x >= 1])
solve!(p)
@assert abs(p.optval - 6.1464) < TOLERANCE

# Test 38
x = Variable(3, 1)
A = -1.0*[0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
c = [3 2 4]
p = maximize(c*x , [quad_form(x, A) >= -1])
solve!(p)
@assert abs(p.optval - 3.7713) < TOLERANCE

# Test 39
x = Variable(3, 1)
A = [2 -3 5; -2 9 -3; 5 -8 3]
b = [-3; 9; 5]
c = [3 2 4]
d = -3
p = minimize(quad_over_lin(A*x + b, c*x + d))
solve!(p)
@assert abs(p.optval - 17.7831) < TOLERANCE

# Test 40
x = Variable(4, 4)
c = ones(16, 1)
reshaped = reshape(x, 16, 1)
a = [1:16]
p = minimize(c' * reshaped, reshaped >= a)
solve!(p)
@assert abs(p.optval - 136) < TOLERANCE

# Test 41
x = Variable(4, 4)
p = minimize(sum(diag(x)), x >= 2)
solve!(p)
@assert abs(p.optval - 8) < TOLERANCE

# Test 42
m = Variable(4, 5)
c = [m[3, 3] == 4, m >= 1]
p = minimize(norm(m, :fro), c)
solve!(p)
@assert abs(m.value[1, 1] - 1) < TOLERANCE

# Test 43
x = Variable()
y = Variable(4)
p = minimize(x + sum(y), [x == 1, y == 3])
solve!(p)
@assert abs(p.optval - 13) < TOLERANCE

# x = Variable(1)
# p = minimize(x, [eye(2) + x >= ones(2, 2)])
# solve!(p)
# p = minimize(y, [X >= 1, y >= X])
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
