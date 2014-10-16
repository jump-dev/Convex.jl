tests = ["test.jl",
         "test_new.jl",
         "test_exp.jl",
         "test_sdp.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end
