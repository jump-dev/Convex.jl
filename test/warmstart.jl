using Convex, SCS
set_default_solver(SCSSolver(verbose=0, max_iters=1000000))

println("small test to compile and see if we throw warnings when warmstarting what can't be warmstarted.")
n = 10
y = rand(n)
x = Variable(n)

lambda = [100]
problem = minimize(sumsquares(y - x) + lambda[1] * sumsquares(x - 10))
@time solve!(problem, warmstart = true)

println("\nnow try first solving and then warmstarting")
n = 1000
y = rand(n)
x = Variable(n)

lambda = [100]
problem = minimize(sumsquares(y - x) + lambda[1] * sumsquares(x - 10))
@time solve!(problem)
@show problem.optval
lambda[1] = 105
println("this run should be faster than the last and have slightly higher optimal value")
@time solve!(problem, warmstart=true)
@show problem.optval

println("\nnow try with scalar variables")
lambda = 100
problem = minimize(sumsquares(y - x) + lambda * sumsquares(x - 10))
@time solve!(problem)
@show problem.optval
lambda = 105
println("this run should be faster than the last and have slightly higher optimal value")
@time solve!(problem, warmstart=true)
@show problem.optval


