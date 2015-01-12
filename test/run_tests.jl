using Convex

tests = ["test_utilities.jl",
         "test_affine.jl",
         "test_lp.jl",
         "test_socp.jl"]
tests_sdp = ["test_sdp.jl"]
tests_exp = ["test_exp.jl"]
tests_int = ["test_int.jl"]
tests_exp_and_sdp = ["test_exp_and_sdp.jl"]

println("Running tests:")

# The following syntax can be used to solve it using other solvers
# using Gurobi
# set_default_solver(GurobiSolver())

for curtest in tests
    info(" Test: $(curtest)")
    include(curtest)
end

if can_solve_sdp(get_default_solver())
    for curtest in tests_sdp
        info(" Test: $(curtest)")
        include(curtest)
    end
end

if can_solve_exp(get_default_solver())
    for curtest in tests_exp
        info(" Test: $(curtest)")
        include(curtest)
    end
end

if can_solve_sdp(get_default_solver()) && can_solve_exp(get_default_solver())
    for curtest in tests_exp_and_sdp
        info(" Test: $(curtest)")
        include(curtest)
    end
end

if can_solve_mip(get_default_solver())
	for curtest in tests_int
    info(" Test: $(curtest)")
    include(curtest)
	end
end
