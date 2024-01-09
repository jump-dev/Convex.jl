using Test

# Seed random number stream to improve test reliability
import Random
Random.seed!(2)

@testset "Convex" begin
    include("test_utilities.jl")
    include("test_problem_depot.jl")
    include("MOI_wrapper.jl")
end
