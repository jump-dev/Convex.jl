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

if isdir(Pkg.dir("SCS"))
	for curtest in tests_scs
    info(" Test: $(curtest)")
    include(curtest)
	end
end

if isdir(Pkg.dir("GLPK")) && isdir(Pkg.dir("GLPKMathProgInterface"))
	for curtest in tests_glpk
    info(" Test: $(curtest)")
    include(curtest)
	end
end
