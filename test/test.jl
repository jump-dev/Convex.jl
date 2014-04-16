@everywhere require("src/CVX.jl")
@everywhere using CVX

n = 2
m = 2
p = 0
l = 2
ncones = 0
q = convert(Ptr{Int64}, C_NULL)
Gpr = [-1.0, -2.0, -2.0, -1.0]
Gjc = [0, 2, 4]
Gir = [0, 1, 0, 1]
Apr = convert(Ptr{Float64}, C_NULL)
Ajc = convert(Ptr{Int64}, C_NULL)
Air = convert(Ptr{Int64}, C_NULL)
c = [1.0, 1.0]
h = [-2.0, -2.0]
b = convert(Ptr{Float64}, C_NULL)
pwork = ccall((:ECOS_setup, "../ecos/ecos.so"), Ptr{Void}, 
	          (Int64, Int64, Int64, Int64, Int64, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, 
	           Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
	           n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)
ccall((:ECOS_solve, "../ecos/ecos.so"), Int64, (Ptr{Void},), pwork)
dptr = convert(Ptr{Ptr{Float64}}, pwork)
ptr = unsafe_load(dptr, 12)
x1 = unsafe_load(ptr, 1)
x2 = unsafe_load(ptr, 2)


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
