# Tests of the write_to_file(problem, filename) interface

# Simple quadratic program
# Can be written in any of the formats
let
    moi_file_formats = ["cbf", "lp", "mof", "mps", "nl", "rew", "sdpa"]
    filename, _io = Base.Filesystem.mktemp()

    m, n = 5, 10
    x = Variable(n)
    A = randn(m, n)
    b = randn(m)
    problem = minimize(sumsquares(A * x - b), [0 <= x, x <= 1])

    # Haven't solved the problem yet
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    # Make sure files are written without error
    # Can't check correctness since we don't provide a read_to_file function
    solve!(problem, SCS.Optimizer)
    for ext in moi_file_formats
        @test write_to_file(problem, filename * ext) |> isnothing
    end
end

# Eric Hanson's fancy problem from issue #395
# https://github.com/jump-dev/Convex.jl/issues/395
let
    filename, _io = Base.Filesystem.mktemp()

    x = HermitianSemidefinite(2)
    problem = minimize(real(tr(x)), tr(x * [1.0 im; -im 0]) == 0, x[1, 1] == 1)

    # Haven't solved the problem yet
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    # Make sure file is written without error
    # Can't check correctness since we don't provide a read_to_file function
    solve!(problem, SCS.Optimizer)
    @test write_to_file(problem, filename * ".sdpa") |> isnothing
end
