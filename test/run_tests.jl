using Convex

tests = ["test.jl",
         "test2.jl"]
tests_scs = ["test_exp.jl",
         "test_sdp.jl"]
tests_glpk = ["test_int.jl"]

println("Running tests:")

for curtest in tests
    info(" Test: $(curtest)")
    include(curtest)
end

if can_solve_sdp(get_default_solver())
	for curtest in tests_scs
    info(" Test: $(curtest)")
    include(curtest)
	end
end

# The following syntax can be used to solve it using other solvers
# using Gurobi
# set_default_solver(GurobiSolver())

if can_solve_mip(get_default_solver())
	for curtest in tests_glpk
    info(" Test: $(curtest)")
    include(curtest)
	end
end
