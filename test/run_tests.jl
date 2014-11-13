tests = ["test.jl",
         "test2.jl"]
tests_scs = ["test_exp.jl",
         "test_sdp.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end

if isdir(Pkg.dir("SCS"))
	for curtest in tests_scs
    println(" Test: $(curtest)")
    include(curtest)
	end
end
