# Tests of the write_to_file(problem, filename) interface

# Simple quadratic program
@testset "Simple quadratic program (Float32 eltype)" begin
    filename, _io = Base.Filesystem.mktemp()

    m, n = 5, 10
    x = Variable(n)
    A = randn(Float32, m, n)
    b = randn(Float32, m)
    problem = minimize(sumsquares(A * x - b), [0 <= x, x <= 1])

    # Haven't solved the problem yet, should get an ArgumentError
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    solve!(problem, SCS.Optimizer)

    # Format lists pulled from source code of MOI.FileFormats.Model
    # These are the formats that are compatible with this problem
    format_ok = [".cbf", ".mof.json", ".dat-s", ".sdpa"]
    for ext in format_ok
        @test write_to_file(problem, filename * ext) |> isnothing
    end

    # These formats are not compatible; the specific exception varies
    format_bad = [".lp", ".mps", ".rew", ".nl"]
    for ext in format_bad
        @test_throws Exception write_to_file(problem, filename * ext)
    end
end

# Simple quadratic program
@testset "Simple quadratic program (BigFloat eltype)" begin
    filename, _io = Base.Filesystem.mktemp()

    m, n = 5, 10
    x = Variable(n)
    A = randn(BigFloat, m, n)
    b = randn(BigFloat, m)
    problem = minimize(sumsquares(A * x - b), [0 <= x, x <= 1])

    # Haven't solved the problem yet, should get an ArgumentError
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    solve!(problem, SCS.Optimizer)

    # Format lists pulled from source code of MOI.FileFormats.Model
    # These are the formats that are compatible with this problem
    format_ok = [".cbf", ".mof.json", ".dat-s", ".sdpa"]
    for ext in format_ok
        @test write_to_file(problem, filename * ext) |> isnothing
    end

    # These formats are not compatible; the specific exception varies
    format_bad = [".lp", ".mps", ".rew", ".nl"]
    for ext in format_bad
        @test_throws Exception write_to_file(problem, filename * ext)
    end
end

# Eric Hanson's fancy problem from issue #395
# https://github.com/jump-dev/Convex.jl/issues/395
@testset "Fancy SDP (Float64 eltype)" begin
    filename, _io = Base.Filesystem.mktemp()

    x = HermitianSemidefinite(2)
    problem = minimize(real(tr(x)), tr(x * [1.0 im; -im 0]) == 0, x[1, 1] == 1)

    # Haven't solved the problem yet, should get an ArgumentError
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    solve!(problem, SCS.Optimizer)

    # Format lists pulled from source code of MOI.FileFormats.Model
    # These are the formats that are compatible with this problem
    format_ok = [".cbf", ".mof.json", ".dat-s", ".sdpa"]
    for ext in format_ok
        @test write_to_file(problem, filename * ext) |> isnothing
    end

    # These formats are not compatible; the specific exception varies
    format_bad = [".lp", ".mps", ".rew", ".nl"]
    for ext in format_bad
        @test_throws Exception write_to_file(problem, filename * ext)
    end
end
