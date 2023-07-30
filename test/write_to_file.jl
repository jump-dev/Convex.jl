# Tests of the write_to_file(problem, filename) interface

# Simple quadratic program
@testset "Simple quadratic program (Float64 eltype)" begin
    filename, _io = Base.Filesystem.mktemp()

    m, n = 5, 10
    x = Variable(n)
    A = rand(m, n)
    b = rand(m)
    problem::Problem{Float64} =
        minimize(sumsquares(A * x - b), [0 <= x, x <= 1])

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

# Simple quadratic program, different eltype
@testset "Simple quadratic program (Float32 eltype)" begin
    filename, _io = Base.Filesystem.mktemp()

    m, n = 5, 10
    x = Variable(n)
    A = rand(Float32, m, n)
    b = rand(Float32, m)
    problem::Problem{Float32} = minimize(
        sumsquares(A * x - b),
        [0 <= x, x <= 1],
        numeric_type = Float32,
    )

    # Haven't solved the problem yet, should get an ArgumentError
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    # Throws an error because SCS doesn't support this eltype
    # (Hypatia.jl gives a similar error)
    @test_broken solve!(problem, SCS.Optimizer)
    return

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

# Same quadratic program, different eltype
@testset "Simple quadratic program (BigFloat eltype)" begin
    filename, _io = Base.Filesystem.mktemp()

    m, n = 5, 10
    x = Variable(n)
    A = rand(BigFloat, m, n)
    b = rand(BigFloat, m)
    problem::Problem{BigFloat} = minimize(
        sumsquares(A * x - b),
        [0 <= x, x <= 1],
        numeric_type = BigFloat,
    )

    # Haven't solved the problem yet, should get an ArgumentError
    @test_throws ArgumentError write_to_file(problem, filename * ".sdpa")

    # Throws an error because SCS doesn't support this eltype
    # (Hypatia.jl gives a similar error)
    @test_broken solve!(problem, SCS.Optimizer)
    return

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
    problem::Problem{Float64} =
        minimize(real(tr(x)), tr(x * [1.0 im; -im 0]) == 0, x[1, 1] == 1)

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
