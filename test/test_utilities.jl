module TestUtilities

using Convex
using Test

import AbstractTrees
import LinearAlgebra
import MathOptInterface as MOI
import SCS
import SparseArrays

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

# It's not super easy to capture the output
# I ended up using this pattern from Suppressor:
# https://github.com/JuliaIO/Suppressor.jl/blob/b4ff08f0fe795a2ce9e592734a758c9e6d8e2bc4/src/Suppressor.jl#L124-L152
function solve_and_return_output(problem, solver; kwargs...)
    original_stdout = stdout
    rd, wr = redirect_stdout()
    out_task = @async read(rd, String)
    try
        solve!(problem, solver; kwargs...)
    finally
        Base.Libc.flush_cstdio() #  https://github.com/JuliaLang/julia/issues/31236
        redirect_stdout(original_stdout)
        close(wr)
    end
    return fetch(out_task)
end

function test_solve!_does_not_return_anything()
    x = Variable()
    p = satisfy(x >= 0)
    output = solve!(
        p,
        MOI.OptimizerWithAttributes(
            SCS.Optimizer,
            "verbose" => 0,
            "eps_abs" => 1e-6,
        ),
    )
    @test output === nothing
    return
end

function test_silent_solver_works()
    x = Variable()
    p = satisfy(x >= 0)
    output_non_silent = solve_and_return_output(
        p,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-6),
    )
    @test output_non_silent != ""
    output_silent = solve_and_return_output(
        p,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-6),
        silent_solver = true,
    )
    @test output_silent == ""
    return
end

function test_solve!_can_take_an_optimizer_directly()
    x = Variable()
    p = satisfy(x >= 0)
    output = solve!(
        p,
        MOI.OptimizerWithAttributes(
            SCS.Optimizer,
            "verbose" => 0,
            "eps_abs" => 1e-6,
        ),
    )
    @test output === nothing
    return
end

function test_complex_objective_function_errors()
    x = Variable()
    @test_throws ErrorException minimize(x + im * x)
    return
end

function test_satisfy_constant_objective()
    x = Variable()
    p = satisfy(x == 0, x == 1)
    @test isnothing(p.objective)
    p = satisfy(Constraint[])
    @test isnothing(p.objective)
    return
end

function test_invalid_head()
    p = Problem(:invalid, nothing, Constraint[])
    err = ErrorException("Unknown type of problem $(p.head)")
    @test_throws err Convex.objective_vexity(p)
    return
end

function test_optval_is_nothing_before_solve!()
    x = Variable()
    p = minimize(x, x >= 0)
    @test p.optval === nothing
    solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    @test p.optval ≈ 0.0 atol = 1e-3
    @test Convex.termination_status(p) == MOI.OPTIMAL
    @test Convex.objective_value(p) ≈ 0.0 atol = 1e-3
    return
end

function test_default_problem_type_is_Float64()
    x = Variable()
    p = minimize(x, x >= 0)
    @test p isa Convex.Problem{Float64}
    return
end

function test_set_value!_doesnt_convert_to_Float64()
    x = Variable()
    set_value!(x, big"1.0")
    @test evaluate(x) isa BigFloat
    x = Variable(2)
    set_value!(x, big.([1.0, 2.0]))
    @test evaluate(x) isa Vector{BigFloat}
    x = Variable(2, 2)
    set_value!(x, big.([1.0 2.0; 3.0 4.0]))
    @test evaluate(x) isa Matrix{BigFloat}
    return
end

function test_show()
    x = Variable()
    @test sprint(show, x) == """
    Variable
    size: (1, 1)
    sign: real
    vexity: affine
    $(Convex.show_id(x))"""
    fix!(x, 1.0)
    @test sprint(show, x) == """
    Variable
    size: (1, 1)
    sign: real
    vexity: constant
    $(Convex.show_id(x))
    value: 1.0"""

    @test sprint(show, 2 * x) == """
        * (constant; real)
        ├─ $(reshape([2], 1, 1))
        └─ real variable (fixed) ($(Convex.show_id(x)))"""

    free!(x)
    p = maximize(log(x), x >= 1, x <= 3)
    @test monotonicity(p) == (Convex.Nonincreasing(),)
    @test sign(p) == NoSign()
    @test curvature(p) == Convex.ConvexVexity()

    @test sprint(show, p) == """
    maximize
    └─ log (concave; real)
       └─ real variable ($(Convex.show_id(x)))
    subject to
    ├─ ≥ constraint (affine)
    │  ├─ real variable ($(Convex.show_id(x)))
    │  └─ $(reshape([1], 1, 1))
    └─ ≤ constraint (affine)
       ├─ real variable ($(Convex.show_id(x)))
       └─ $(reshape([3], 1, 1))

    status: `solve!` not called yet"""

    x = ComplexVariable(2, 3)
    @test sprint(show, x) == """
    Variable
    size: (2, 3)
    sign: complex
    vexity: affine
    $(Convex.show_id(x))"""

    # test `MAXDEPTH`
    # We construct a binary tree of depth >= 3
    # to make sure it gets truncated appropriately.
    x = Variable(2)
    y = Variable(2)
    level3 = hcat(x, y)
    level2 = hcat(level3, level3)
    root = hcat(level2, level2)
    p = minimize(sum(x), root == root)
    @test curvature(p) == Convex.ConstVexity()
    @test sprint(show, p) == """
    minimize
    └─ sum (affine; real)
       └─ 2-element real variable ($(Convex.show_id(x)))
    subject to
    └─ == constraint (affine)
       ├─ hcat (affine; real)
       │  ├─ hcat (affine; real)
       │  │  ├─ …
       │  │  └─ …
       │  └─ hcat (affine; real)
       │     ├─ …
       │     └─ …
       └─ hcat (affine; real)
          ├─ hcat (affine; real)
          │  ├─ …
          │  └─ …
          └─ hcat (affine; real)
             ├─ …
             └─ …

    status: `solve!` not called yet"""

    # test `MAXWIDTH`
    x = Variable()
    p = satisfy([x == i for i in 1:100])
    err = ErrorException("Satisfiability problem cannot be used as subproblem")
    @test_throws err sign(p)
    old_maxwidth = Convex.MAXWIDTH[]
    Convex.MAXWIDTH[] = 2
    @test sprint(show, p) == """
        satisfy
        └─ nothing
        subject to
        ├─ == constraint (affine)
        │  ├─ real variable ($(Convex.show_id(x)))
        │  └─ $(reshape([1], 1, 1))
        ├─ == constraint (affine)
        │  ├─ real variable ($(Convex.show_id(x)))
        │  └─ $(reshape([2], 1, 1))
        ⋮

        status: `solve!` not called yet"""
    Convex.MAXWIDTH[] = old_maxwidth

    # solved problem
    x = Variable()
    p = satisfy(x >= 0)
    output = solve!(
        p,
        MOI.OptimizerWithAttributes(
            SCS.Optimizer,
            "verbose" => 0,
            "eps_abs" => 1e-6,
        ),
    )
    @test sprint(show, p) == """
            satisfy
            └─ nothing
            subject to
            └─ ≥ constraint (affine)
               ├─ real variable ($(Convex.show_id(x)))
               └─ $(reshape([0], 1, 1))

            termination status: OPTIMAL
            primal status: FEASIBLE_POINT
            dual status: FEASIBLE_POINT"""

    # test small `MAXDIGITS`
    x = Variable()
    old_maxdigits = Convex.MAXDIGITS[]
    Convex.MAXDIGITS[] = 2
    @test length(Convex.show_id(x)) == length("id: ") + 5
    Convex.MAXDIGITS[] = old_maxdigits

    # test large `MAXDIGITS`
    x = Variable()
    old_maxdigits = Convex.MAXDIGITS[]
    Convex.MAXDIGITS[] = 100
    @test length(Convex.show_id(x)) ==
          length("id: ") + length(string(x.id_hash))
    Convex.MAXDIGITS[] = old_maxdigits
    return
end

function test_vartype_and_set_vartype()
    for x in (Variable(), Variable(1))
        @test vartype(x) == ContVar

        vartype!(x, BinVar)
        @test vartype(x) == BinVar
        @test x.vartype == BinVar

        vartype!(x, IntVar)
        @test vartype(x) == IntVar
        @test x.vartype == IntVar

        vartype!(x, ContVar)
        @test vartype(x) == ContVar
        @test x.vartype == ContVar
    end
    return
end

function test_Constructors()
    # Constructors with sign
    for sgn in (Positive(), NoSign())
        for x in [   # tuple size
            Variable((2, 2), sgn),
            Variable((2, 2), sgn, BinVar),
            Variable((2, 2), sgn, :Bin),
            # individual size
            Variable(2, 2, sgn),
            Variable(2, 2, sgn, BinVar),
            Variable(2, 2, sgn, :Bin),
            # single dimension
            Variable(2, sgn),
            Variable(2, sgn, BinVar),
            Variable(2, sgn, :Bin),
            # no dimension
            Variable(sgn),
            Variable(sgn, BinVar),
            Variable(sgn, :Bin),
        ]
            @test x isa Variable
            @test x isa Convex.AbstractVariable
            @test sign(x) == sgn
            @test x.sign == sgn
        end
    end

    # constructors without sign
    for x in [   # tuple size
        Variable((2, 2)),
        Variable((2, 2), BinVar),
        Variable((2, 2), :Bin),
        # individual size
        Variable(2, 2),
        Variable(2, 2, BinVar),
        Variable(2, 2, :Bin),
        # single dimension
        Variable(2),
        Variable(2, BinVar),
        Variable(2, :Bin),
        # no dimension
        Variable(),
        Variable(BinVar),
        Variable(:Bin),
    ]
        @test x isa Variable
        @test x isa Convex.AbstractVariable
        @test sign(x) == NoSign()
        @test x.sign == NoSign()

        Convex.sign!(x, Positive())
        @test sign(x) == Positive()
        @test x.sign == Positive()
    end

    # ComplexVariable
    for x in [   # tuple size
        ComplexVariable((2, 2)),
        Variable((2, 2), ComplexSign()),
        ComplexVariable((2, 2), :Semidefinite),
        # individual size
        ComplexVariable(2, 2),
        Variable(2, 2, ComplexSign()),
        ComplexVariable(2, 2, :Semidefinite),
        # single dimension
        ComplexVariable(2),
        Variable(2, ComplexSign()),
        # no dimension
        ComplexVariable(),
        Variable(ComplexSign()),
    ]
        @test x isa ComplexVariable
        @test x isa Convex.AbstractVariable
        @test sign(x) == ComplexSign()
    end

    for vt in (BinVar, IntVar),
        V in (ComplexVariable, Semidefinite, HermitianSemidefinite)

        @test_throws Any V(2; vartype = vt)
    end

    for vt in (:Bin, :Int),
        V in (Semidefinite, HermitianSemidefinite, ComplexVariable)

        @test_throws Any V(2, vt)
    end

    # Semidefinite
    for x in [
        Variable((2, 2), :Semidefinite),
        Variable(2, 2, :Semidefinite),
        ComplexVariable((2, 2), :Semidefinite),
        ComplexVariable(2, 2, :Semidefinite),
        HermitianSemidefinite((2, 2)),
        HermitianSemidefinite(2, 2),
        HermitianSemidefinite(2),
        Semidefinite((2, 2)),
        Semidefinite(2, 2),
        Semidefinite(2),
    ]
        @test length(get_constraints(x)) == 1
        @test get_constraints(x)[] isa Convex.PositiveSemidefiniteConeConstraint
    end

    @test_throws ErrorException HermitianSemidefinite(2, 3)
    @test_throws ErrorException Semidefinite(2, 3)
    return
end

function test_length_and_size()
    x = Variable(2, 3)
    @test length(x) == 6
    @test size(x) == (2, 3)
    @test size(x, 1) == 2
    @test size(x, 2) == 3

    x = Variable(3)
    @test length(x) == 3
    @test size(x) == (3, 1)

    x = Variable()
    @test length(x) == 1
    @test size(x) == (1, 1)
    return
end

function test_lastindex_and_axes()
    x = Variable(2, 3)
    @test axes(x) == (Base.OneTo(2), Base.OneTo(3))
    @test axes(x, 1) == Base.OneTo(2)
    @test lastindex(x) == 6
    @test lastindex(x, 2) == 3
    y = x[:, end]
    @test y isa Convex.AbstractExpr
    @test size(y) == (2, 1)
    return
end

function test_Cartesian_index()
    x = Variable(3, 2)
    set_value!(x, rand(3, 2))
    context = Convex.Context{Float64}(() -> MOI.Utilities.Model{Float64}())

    for ind in CartesianIndices(zeros(3, 2))
        L = Convex.conic_form!(context, x[ind])
        R = Convex.conic_form!(context, x[ind[1], ind[2]])
        @test L == R
        @test evaluate(x[ind]) == evaluate(x[ind[1], ind[2]])
    end
    y = [1.0 2 3; 4 5 6] * x
    for ind in CartesianIndices(zeros(2, 2))
        L = Convex.conic_form!(context, y[ind])
        R = Convex.conic_form!(context, y[ind[1], ind[2]])
        @test L == R
        @test evaluate(y[ind]) == evaluate(y[ind[1], ind[2]])
    end
    return
end

function test_parametric_constants()
    z = constant([1.0 0.0im; 0.0 1.0])
    @test z isa Convex.ComplexConstant{Float64}
    # Helper functions
    @test Convex.ispos(1)
    @test Convex.ispos(0)
    @test !Convex.ispos(-1)
    @test Convex.ispos([0, 1, 0])
    @test !Convex.ispos([0, -1, 0])
    @test Convex.isneg(-1)
    @test Convex.isneg(0)
    @test !Convex.isneg(1)
    @test Convex.isneg([0, -1, 0])
    @test !Convex.isneg([0, 1, 0])
    @test Convex._size(3) == (1, 1)
    @test Convex._sign(3) == Positive()
    @test Convex._size([-1, 1, 1]) == (3, 1)
    @test Convex._sign([-1, 1, 1]) == NoSign()
    @test Convex._sign([-1, -1, -1]) == Negative()
    @test Convex._size([0 0; 0 0]) == (2, 2)
    @test Convex._sign([0 0; 0 0]) == Positive()
    @test Convex._size(0 + 1im) == (1, 1)
    @test Convex._sign(0 + 1im) == ComplexSign()
    return
end

function test_issue_341_evaluate_for_constants()
    A = rand(4, 4)
    @test evaluate(constant(A)) ≈ copy(A)
    @test constant(A).size == (4, 4)
    b = rand(4)
    @test evaluate(constant(b)) ≈ copy(b)
    @test constant(b).size == (4, 1)
    c = 1.0
    @test evaluate(constant(c)) ≈ c
    @test constant(c).size == (1, 1)

    @test evaluate(sumlargesteigs(Variable(4, 4), 0)) == 0
    @test evaluate(sumlargest(Variable(4), 0)) == 0
    @test evaluate(sumsmallest(Variable(4), 0)) == 0
    return
end

function test_Base_vect()
    # Issue #223: ensure we can make vectors of variables
    @test size([Variable(2), Variable(3, 4)]) == (2,)
    return
end

function test_Iteration()
    x = Variable(2, 3)
    s = sum([xi for xi in x])
    set_value!(x, [1 2 3; 4 5 6])
    # evaluate(s) == [21] (which might be wrong? expected 21)
    # but [21][1] === 21[1] === 21
    # so this should pass even after "fixing" that
    @test evaluate(s)[1] == 21

    x = Variable(4)
    @test [xi.inds for xi in x] == [1:1, 2:2, 3:3, 4:4]

    x = Variable(0)
    @test [xi for xi in x] == []
    @test iterate(x) == nothing
    return
end
# returns [21]; not sure why
# context("iteration") do
#     x = Variable(2,3)
#     s = sum([xi for xi in x])
#     x.value = [1 2 3; 4 5 6]
#     @fact evaluate(s) --> 21
# end

function test_DCP_warnings()
    x = Variable()
    y = Variable()
    p = minimize(log(x) + square(y), x >= 0, y >= 0)
    @test_throws DCPViolationError solve!(p, SCS.Optimizer)
    str = sprint(Base.showerror, DCPViolationError())
    @test contains(str, "Expression not DCP compliant")

    p = minimize(sqrt(x), x >= 0, x <= 1)
    @test_throws DCPViolationError solve!(p, SCS.Optimizer)
    return
end

function test_add_constraints!_issue_380()
    x = Variable(3, 3)
    p = minimize(norm(x, 1))
    y = randn(3, 3)
    c = (norm2(x - y) <= 1)
    @test length(p.constraints) == 0
    add_constraint!(p, c)
    @test length(p.constraints) == 1
    empty!(p.constraints)
    add_constraints!(p, c)
    @test length(p.constraints) == 1
    empty!(p.constraints)
    add_constraint!(p, [c])
    @test length(p.constraints) == 1
    empty!(p.constraints)
    c2 = (norm2(x - rand(3, 3)) <= 3)
    add_constraints!(p, [c, c2])
    @test length(p.constraints) == 2
    return
end

function test_diagm_issue_401()
    x = Variable(3)
    @test diagm(x) isa Convex.AbstractExpr
    return
end

function test_is_psd()
    _test_is_psd(Float64)
    _test_is_psd(Float32)
    _test_is_psd(Int)
    _test_is_psd(ComplexF64)
    _test_is_psd(BigFloat)
    _test_is_psd(Rational{BigInt})
    _test_is_psd(Complex{Rational{BigInt}})
    return
end

function _test_is_psd(T)
    A = zeros(T, 3, 3)
    A[1, 1] = one(T)
    @test Convex._is_psd(A)
    @test Convex._is_psd(SparseArrays.sparse(A))
    B = A .- one(T) / T(5000)
    @test !Convex._is_psd(B)
    @test !Convex._is_psd(SparseArrays.sparse(B))

    # See https://github.com/jump-dev/Convex.jl/issues/452 for details
    C = [
        70.12718378756115 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        -70.12718378756115 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 277.05985725905305 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 70.12718378756115 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 347.1870410466142 -103.46633673574595 0.0 -70.12718378756115 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 103.46633673574595 -103.46633673574595 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -70.12718378756115 0.0 -103.46633673574595 277.05985725905305 -103.46633673574595
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -103.46633673574595 103.46633673574595
    ]
    @test Convex._is_psd(C)
    @test Convex._is_psd(SparseArrays.sparse(C))
    return
end

function test_assume_psd_option_for_quadform()
    A = [-1 0; 0 1] # neither PSD nor negative semidefinite
    x = Variable(2)

    @test_throws ErrorException quadform(x, A) # default
    @test_throws ErrorException quadform(x, A; assume_psd = false)
    @test quadform(x, A; assume_psd = true) isa Convex.AbstractExpr
    return
end

function test_logsumexp_stability()
    v = Convex.constant([1000, 1000, 1000])
    @test Convex.evaluate(Convex.logsumexp(v)) ≈ 1001.098612
    return
end

function test_conv_issue_364()
    n = 3
    m = 11
    h = rand(m)
    x = rand(n)
    hvar = Variable(m)
    hvar.value = h
    function _conv(h, x)
        m = length(h)
        n = length(x)
        zero_pad_x(i) = 1 <= i <= n ? x[i] : 0
        return [sum(h[j] * zero_pad_x(i - j + 1) for j in 1:m) for i in 1:m+n-1]
    end
    @test evaluate(conv(hvar, x)) ≈ _conv(h, x)
    return
end

function test_conj_issue_416()
    A = [1 1im; -1im 1]
    X = ComplexVariable(2, 2)
    p = minimize(real(tr(conj(X))), [X == A])
    solve!(
        p,
        MOI.OptimizerWithAttributes(
            SCS.Optimizer,
            "verbose" => 1,
            "eps_abs" => 1e-6,
        ),
    )
    @test evaluate(X) ≈ A atol = 1e-3
    return
end

function test_logisticloss_issue_458()
    x = Variable()
    expr = logisticloss(x)
    @test expr isa Convex.AbstractExpr
    set_value!(x, 1.5)
    @test evaluate(expr) ≈ log(1 + exp(1.5))
    return
end

function test_dot_issue_508()
    x = [1.0 + 1.0im]
    y = [-1.0im]
    @test dot(x, y) ≈ evaluate(dot(constant(x), y))
    return
end

function _test_sparse_tape(T)
    d_in = 5
    variables = MOI.VariableIndex.(1:d_in)
    input = rand(T, d_in)
    A = SparseArrays.sprand(T, d_in, d_in, 0.1)
    b = SparseArrays.sprand(T, d_in, 0.8)
    A_init = copy(A)
    b_init = copy(b)
    op = Convex.SparseAffineOperation(A, b)
    tape = Convex.SparseTape(op, variables)
    collapsed_tape = Convex.SparseAffineOperation(tape)
    @test collapsed_tape.matrix * input + collapsed_tape.vector ≈ A * input + b
    op2 = Convex.SparseAffineOperation(
        SparseArrays.sparse(one(T) * LinearAlgebra.I, d_in, d_in),
        -b,
    )
    tape = Convex.add_operation(tape, op2)
    collapsed_tape2 = Convex.SparseAffineOperation(tape)
    @test collapsed_tape2.matrix * input + collapsed_tape2.vector ≈ A * input
    op3 = Convex.SparseAffineOperation(ones(T, 1, d_in), [zero(T)])
    tape = Convex.add_operation(tape, op3)
    collapsed_tape3 = Convex.SparseAffineOperation(tape)
    @test collapsed_tape3.matrix * input + collapsed_tape3.vector ≈
          [sum(A * input)]
    @test A_init ≈ A
    @test b_init ≈ b
    return
end

test_sparse_tape_Float64() = _test_sparse_tape(Float64)

test_sparse_tape_Float32() = _test_sparse_tape(Float32)

test_sparse_tape_BigFloat() = _test_sparse_tape(BigFloat)

module DictVectors

using Convex

# To make sure `Convex` isn't using field access on `AbstractVariable`'s
# we'll use a global dictionary to store information about each instance
# our of mock variable type, `DictVector`.
const global_cache = Dict{UInt64,Any}()

mutable struct DictVector{T} <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int,Int}
    function DictVector{T}(d) where {T}
        this = new(:DictVector, 0, (d, 1))
        this.id_hash = objectid(this)
        global_cache[this.id_hash] = Dict(
            :value => nothing,
            :sign => T <: Complex ? ComplexSign() : NoSign(),
            :vartype => ContVar,
            :constraints => Constraint[],
            :vexity => Convex.AffineVexity(),
        )
        return this
    end
end

Convex.evaluate(x::DictVector) = global_cache[x.id_hash][:value]

function Convex.set_value!(x::DictVector, v::AbstractArray)
    return global_cache[x.id_hash][:value] = v
end

function Convex.set_value!(x::DictVector, v::Number)
    return global_cache[x.id_hash][:value] = v
end

Convex.vexity(x::DictVector) = global_cache[x.id_hash][:vexity]

function Convex.vexity!(x::DictVector, v::Convex.Vexity)
    return global_cache[x.id_hash][:vexity] = v
end

Convex.sign(x::DictVector) = global_cache[x.id_hash][:sign]

Convex.sign!(x::DictVector, s::Convex.Sign) = global_cache[x.id_hash][:sign] = s

Convex.vartype(x::DictVector) = global_cache[x.id_hash][:vartype]

function Convex.vartype!(x::DictVector, s::Convex.VarType)
    return global_cache[x.id_hash][:vartype] = s
end

Convex.get_constraints(x::DictVector) = global_cache[x.id_hash][:constraints]

function Convex.add_constraint!(x::DictVector, s::Convex.Constraint)
    return push!(global_cache[x.id_hash][:constraints], s)
end

end  # DictVector

function test_DictVectors()
    # Let us solve a basic problem from `test_affine.jl`
    x = DictVectors.DictVector{BigFloat}(1)
    y = DictVectors.DictVector{BigFloat}(1)
    p = minimize(x + y, [x >= 3, y >= 2])
    @test vexity(p) == Convex.AffineVexity()
    solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    @test p.optval ≈ 5 atol = 1e-3
    @test evaluate(x + y) ≈ 5 atol = 1e-3

    add_constraint!(x, x >= 4)
    solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    @test p.optval ≈ 6 atol = 1e-3
    @test evaluate(x + y) ≈ 6 atol = 1e-3
    @test length(get_constraints(x)) == 1
    return
end

module DensityMatricies

using Convex

mutable struct DensityMatrix{T} <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int,Int}
    value::Union{Convex.Value,Nothing}
    vexity::Convex.Vexity
    function DensityMatrix(d)
        this = new{ComplexF64}(
            :DensityMatrix,
            0,
            (d, d),
            nothing,
            Convex.AffineVexity(),
        )
        this.id_hash = objectid(this)
        return this
    end
end

Convex.get_constraints(ρ::DensityMatrix) = [ρ ⪰ 0, tr(ρ) == 1]

Convex.sign(::DensityMatrix) = Convex.ComplexSign()

Convex.vartype(::DensityMatrix) = Convex.ContVar

end  # DensityMatricies

function test_DensityMatrix()
    X = randn(4, 4) + im * rand(4, 4)
    X = X + X'
    # `X` is Hermitian and non-degenerate (with probability 1)
    # Let us calculate the projection onto the eigenspace associated
    # to the maximum eigenvalue
    e_vals, e_vecs = LinearAlgebra.eigen(LinearAlgebra.Hermitian(X))
    e_val, idx = findmax(e_vals)
    e_vec = e_vecs[:, idx]
    proj = e_vec * e_vec'
    # found it! Now let us do it again via an SDP
    ρ = DensityMatricies.DensityMatrix(4)
    prob = maximize(real(tr(ρ * X)))
    solve!(prob, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    @test prob.optval ≈ e_val atol = 1e-3
    @test evaluate(ρ) ≈ proj atol = 1e-3
    ρ2 = real(ρ) + im * imag(ρ)
    prob = maximize(real(tr(ρ2 * X)))
    solve!(prob, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    @test prob.optval ≈ e_val atol = 1e-3
    @test evaluate(ρ) ≈ proj atol = 1e-3
    @test evaluate(ρ) ≈ evaluate(ρ2) atol = 1e-3
    return
end

module ProbabilityVectors

using Convex

mutable struct ProbabilityVector <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int,Int}
    value::Union{Convex.Value,Nothing}
    vexity::Convex.Vexity
    function ProbabilityVector(d)
        this =
            new(:ProbabilityVector, 0, (d, 1), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        return this
    end
end

Convex.get_constraints(p::ProbabilityVector) = [sum(p) == 1]

Convex.sign(::ProbabilityVector) = Convex.Positive()

Convex.vartype(::ProbabilityVector) = Convex.ContVar

(p::ProbabilityVector)(x) = dot(p, x)

end  # ProbabilityVectors

function test_ProbabilityVectors()
    p = ProbabilityVectors.ProbabilityVector(3)
    x = [1.0, 2.0, 3.0]
    @test p(x) isa Convex.AbstractExpr
    @test sign(p) == Positive()
    prob = minimize(p(x))
    solve!(prob, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    @test prob.optval ≈ 1.0 atol = 1e-3
    @test evaluate(p(x)) ≈ 1.0 atol = 1e-3
    @test evaluate(p) ≈ [1.0, 0.0, 0.0] atol = 1e-3
    return
end

function test_write_to_file()
    x = Variable(3)
    p = minimize(logsumexp(x))
    dir = mktempdir()
    filename = joinpath(dir, "test.mof.json")
    @test_throws ArgumentError write_to_file(p, filename)
    solve!(p, SCS.Optimizer; silent_solver = true)
    write_to_file(p, filename)
    @test occursin("ExponentialCone", read(filename, String))
    p_int = minimize(logsumexp(x); numeric_type = Int)
    @test_throws MethodError write_to_file(p_int, filename)
    return
end

function test_deprecation_strict_inequality()
    x = Variable()
    @test_logs (:warn,) x < 1
    @test_logs (:warn,) 1 < x
    @test_logs (:warn,) x < x
    @test_logs (:warn,) x > 1
    @test_logs (:warn,) 1 > x
    @test_logs (:warn,) x > x
    @test string(x < 1) == string(x <= 1)
    @test string(x > 1) == string(x >= 1)
    return
end

function test_deprecation_norm()
    x = Variable(2)
    @test_deprecated norm_inf(x)
    @test_deprecated norm_1(x)
    @test_deprecated norm_fro(x)
    return
end

function test_deprecation_in_symbol()
    x = Variable(2, 2)
    @test_logs (:warn,) (x in :SDP)
    @test in(x, :semidefinite) isa Convex.PositiveSemidefiniteConeConstraint
    return
end

function test_dcp_rules()
    vexities = (
        Convex.ConcaveVexity(),
        Convex.ConvexVexity(),
        Convex.ConstVexity(),
        Convex.AffineVexity(),
        Convex.NotDcp(),
    )
    signs = (
        Convex.Positive(),
        Convex.Negative(),
        Convex.NoSign(),
        Convex.ComplexSign(),
    )
    monotonicities = (
        Convex.Nondecreasing(),
        Convex.Nonincreasing(),
        Convex.ConstMonotonicity(),
        Convex.NoMonotonicity(),
    )
    # -(::Vexity)
    for (i, j) in enumerate([2, 1, 3, 4, 5])
        @test -vexities[i] == vexities[j]
    end
    # -(::Monotonicity)
    for (i, j) in enumerate([2, 1, 3, 4])
        @test -monotonicities[i] == monotonicities[j]
    end
    # +(::Vexity, ::Vexity)
    add_rule_vexity = [
        1 5 1 1 5
        5 2 2 2 5
        1 2 3 4 5
        1 2 4 4 5
        5 5 5 5 5
    ]
    for i in 1:size(add_rule_vexity, 1), j in 1:size(add_rule_vexity, 2)
        @test vexities[i] + vexities[j] == vexities[add_rule_vexity[i, j]]
    end
    # +(::Sign)
    for i in 1:4
        @test +(signs[i]) == signs[i]
    end
    # -(::Sign)
    for (i, j) in enumerate([2, 1, 3, 4])
        @test -signs[i] == signs[j]
    end
    # +(::Sign, ::Sign)
    add_rule_sign = [
        1 3 3 4
        3 2 3 4
        3 3 3 4
        4 4 4 4
    ]
    for i in 1:size(add_rule_sign, 1), j in 1:size(add_rule_sign, 2)
        @test signs[i] + signs[j] == signs[add_rule_sign[i, j]]
    end
    # *(::Sign, ::Sign)
    mul_rule_sign = [
        1 2 3 4
        2 1 3 4
        3 3 3 4
        4 4 4 4
    ]
    for i in 1:size(mul_rule_sign, 1), j in 1:size(mul_rule_sign, 2)
        @test signs[i] * signs[j] == signs[mul_rule_sign[i, j]]
    end
    # *(::Sign, ::Monotonicity)
    # *(::Monotonicity, ::Sign)
    mul_rule_sign_mon = [
        1 2 3 4
        2 1 3 4
        4 4 4 4
        4 4 4 4
    ]
    for i in 1:size(mul_rule_sign_mon, 1), j in 1:size(mul_rule_sign_mon, 2)
        @test signs[i] * monotonicities[j] ==
              monotonicities[mul_rule_sign_mon[i, j]]
        @test monotonicities[j] * signs[i] ==
              monotonicities[mul_rule_sign_mon[i, j]]
    end
    # *(::Vexity, ::Monotonicity)
    # *(::Monotonicity, ::Vexity)
    mul_rule_mon_vex = [
        1 2 3 4 5
        2 1 3 4 5
        1 2 3 4 5
        5 5 3 4 5
    ]
    for i in 1:size(mul_rule_mon_vex, 1), j in 1:size(mul_rule_mon_vex, 2)
        @test monotonicities[i] * vexities[j] ==
              vexities[mul_rule_mon_vex[i, j]]
        @test vexities[j] * monotonicities[i] ==
              vexities[mul_rule_mon_vex[i, j]]
    end
    # *(::Vexity, ::Sign)
    # *(::Sign, ::Vexity)
    mul_rule_sign_vex = [
        1 2 3 4 5
        1 2 3 4 5
        1 2 3 4 5
        5 5 3 4 5
    ]
    for i in 1:size(mul_rule_sign_vex, 1), j in 1:size(mul_rule_sign_vex, 2)
        @test signs[i] * vexities[j] == vexities[mul_rule_sign_vex[i, j]]
        @test vexities[j] * signs[i] == vexities[mul_rule_sign_vex[i, j]]
    end
    return
end

function test_problem_maximize()
    x = Variable(1, Positive())
    p = maximize(exp(x), x <= 1)
    @test monotonicity(p) == (Convex.Nonincreasing(),)
    @test curvature(p) == Convex.ConcaveVexity()
    @test sign(p) == Convex.Negative()
    return
end

function test_problem_minimize()
    x = Variable(1, Positive())
    p = minimize(exp(x), x <= 1)
    @test monotonicity(p) == (Convex.Nondecreasing(),)
    @test curvature(p) == Convex.ConvexVexity()
    @test sign(p) == Convex.Positive()
    return
end

function test_problem_satisfy()
    p = satisfy(Variable() >= 0)
    @test_throws(
        ErrorException("Satisfiability problem cannot be used as subproblem"),
        monotonicity(p),
    )
    @test_throws(
        ErrorException("Satisfiability problem cannot be used as subproblem"),
        curvature(p),
    )
    @test_throws(
        ErrorException("Satisfiability problem cannot be used as subproblem"),
        sign(p),
    )
    return
end

function test_expressions()
    x = Variable(2)
    @test hash(x) isa UInt
    @test size(x) == (2, 1)
    @test size(x, 1) == 2
    @test size(x, 2) == 1
    @test size(x, 3) == 1
    @test_throws ErrorException("dimension out of range") size(x, -1)
    @test ndims(x) == 2
    @test lastindex(x) == 2
    @test axes(x) == (1:2, 1:1)
    @test axes(x, 1) == 1:2
    @test axes(x, 2) == 1:1
    @test lastindex(x, 1) == 2
    @test lastindex(x, 2) == 1
    return
end

function test_tree_interface()
    x = Variable()
    objective = exp(x)
    constraints = [x >= 0]
    p = minimize(objective, constraints)
    o, c = AbstractTrees.children(p)
    @test o === objective
    @test c[1] === constraints[1]
    @test sprint(AbstractTrees.printnode, constraints) === "constraints"
    @test sprint(AbstractTrees.printnode, p) === "minimize"
    @test sprint(AbstractTrees.printnode, Constraint[]) === "no constraints"
    return
end

function test_multiple_constraint_dual()
    x = Variable()
    c = x >= 1
    p = minimize(x, [c, c])
    solve!(p, SCS.Optimizer)
    @test isapprox(c.dual, 1.0; atol = 1e-5)
    return
end

function test_fixed_variable_value()
    x = Variable()
    y = Variable()
    fix!(x, 2.0)
    p = minimize(y, x + y >= 1)
    solve!(p, SCS.Optimizer)
    @test isapprox(x.value, 2.0; atol = 1e-5)
    @test isapprox(y.value, -1.0; atol = 1e-5)
    return
end

function test_scalar_fn_constant_objective()
    x = Variable()
    p = minimize(2.1, [x >= 1])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(p.optval, 2.1; atol = 1e-5)
    p = minimize(2.2, x >= 1)
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(p.optval, 2.2; atol = 1e-5)
    p = maximize(2.3, [x >= 1])
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(p.optval, 2.3; atol = 1e-5)
    p = maximize(2.4, x >= 1)
    solve!(p, SCS.Optimizer; silent_solver = true)
    @test isapprox(p.optval, 2.4; atol = 1e-5)
    return
end

function test_scalar_fn_objective_number()
    x = Variable()
    p = minimize(constant(2), [x >= 1])
    solve!(p, SCS.Optimizer)
    @test isapprox(p.optval, 2.0; atol = 1e-5)
    return
end

function test_scalar_fn_objective_variable()
    x = Variable()
    p = minimize(x, [x >= 1])
    solve!(p, SCS.Optimizer)
    @test isapprox(p.optval, 1.0; atol = 1e-5)
    return
end

function test_scalar_fn_objective_affine()
    x = Variable()
    p = minimize(x + 1, [x >= 1])
    solve!(p, SCS.Optimizer)
    @test isapprox(p.optval, 2.0; atol = 1e-5)
    return
end

function test_scalar_fn_objective_square()
    x = Variable()
    p = minimize(square(x - 2), [x >= 1])
    solve!(p, SCS.Optimizer)
    @test isapprox(p.optval, 0.0; atol = 1e-3)
    return
end

function test_variable_errors()
    @test_throws ArgumentError Variable((2, 2), ComplexSign(), BinVar)
    return
end

function test_complex_variable_errors()
    y = @test_logs (:warn,) ComplexVariable(:Semidefinite)
    @test size(y) == (1, 1)
    return
end

function test_set_value_errors()
    x = Variable()
    @test_throws DimensionMismatch set_value!(x, [1.0 2.0; 3.0 4.0])
    @test_throws DimensionMismatch set_value!(x, [1.0, 2.0])
    return
end

function test_set_value_nothing()
    x = Variable()
    @test_throws(
        ErrorException("Value of the variable is yet to be calculated"),
        evaluate(x)
    )
    set_value!(x, 1.0)
    @test x.value == 1.0
    set_value!(x, nothing)
    @test x.value === nothing
    return
end

function test_set_value_complex()
    x = ComplexVariable()
    set_value!(x, 1.0)
    @test x.value == 1.0 + 0.0im
    y = ComplexVariable(2)
    set_value!(y, [1.0, 2.0])
    @test y.value == [1.0 + 0.0im, 2.0 + 0.0im]
    return
end


end  # TestUtilities

TestUtilities.runtests()
