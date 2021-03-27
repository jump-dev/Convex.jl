using Convex: AbstractExpr, ConicObj
using LinearAlgebra

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

@testset "Utilities" begin

    @testset "`solve!` does not return anything" begin
        x = Variable()
        p = satisfy(x >= 0)
        output = solve!(p, () -> SCS.Optimizer(verbose=0, eps=1e-6))
        @test output === nothing
    end

    @testset "`silent_solver` works" begin
        x = Variable()
        p = satisfy(x >= 0)
        output_non_silent = solve_and_return_output(p, () -> SCS.Optimizer(eps=1e-6))
        @test output_non_silent != ""
        output_silent = solve_and_return_output(p, () -> SCS.Optimizer(eps=1e-6), silent_solver=true)
        @test output_silent == ""
    end

    # This might get deprecated later.
    @testset "`solve!` can take an optimizer directly" begin
        x = Variable()
        p = satisfy(x >= 0)
        output = solve!(p, SCS.Optimizer(verbose=0, eps=1e-6))
        @test output === nothing
    end

    @testset "`solve!` with MPB solver errors" begin
        x = Variable()
        p = satisfy(x >= 0)
        @test_throws ErrorException solve!(p, SCSSolver())
    end

    @testset "Complex objective function errors" begin
        x = Variable()
        @test_throws ErrorException minimize(x + im*x)
    end

    @testset "`optval` is nothing before `solve!`" begin
        x = Variable()
        p = minimize(x, x >= 0)
        @test p.optval === nothing
        solve!(p, () -> SCS.Optimizer(verbose=false))
        @test p.optval ≈ 0.0 atol=1e-3
        @test Convex.termination_status(p) == MOI.OPTIMAL
        @test Convex.objective_value(p) ≈ 0.0 atol=1e-3
    end

    @testset "Default problem type is `Float64`" begin
        x = Variable()
        p = minimize(x, x >= 0)
        @test p isa Convex.Problem{Float64}
    end

    @testset "`set_value!` doesn't convert to `Float64`" begin
        x = Variable()
        set_value!(x, big"1.0")
        @test evaluate(x) isa BigFloat

        x = Variable(2)
        set_value!(x, big.([1.0, 2.0]))
        @test evaluate(x) isa Vector{BigFloat}

        x = Variable(2, 2)
        set_value!(x, big.([1.0 2.0; 3.0 4.0]))
        @test evaluate(x) isa Matrix{BigFloat}
    end

    @testset "Show" begin
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

        @test sprint(show, 2*x) == """
        * (constant; real)
        ├─ 2
        └─ real variable (fixed) ($(Convex.show_id(x)))"""
        
        free!(x)
        p = maximize( log(x), x >= 1, x <= 3 )

        @test sprint(show, p) == """
        maximize
        └─ log (concave; real)
           └─ real variable ($(Convex.show_id(x)))
        subject to
        ├─ >= constraint (affine)
        │  ├─ real variable ($(Convex.show_id(x)))
        │  └─ 1
        └─ <= constraint (affine)
           ├─ real variable ($(Convex.show_id(x)))
           └─ 3
        
        status: `solve!` not called yet"""
        
        x = ComplexVariable(2,3)
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
        level3 = hcat(x,y)
        level2 = hcat(level3, level3)
        root = hcat(level2, level2)
        p = minimize(sum(x), root == root)
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
        p = satisfy([ x == i for i = 1:100])
        old_maxwidth = Convex.MAXWIDTH[]
        Convex.MAXWIDTH[] = 2
        @test sprint(show, p) == """
            minimize
            └─ 0
            subject to
            ├─ == constraint (affine)
            │  ├─ real variable ($(Convex.show_id(x)))
            │  └─ 1
            ├─ == constraint (affine)
            │  ├─ real variable ($(Convex.show_id(x)))
            │  └─ 2
            ⋮

            status: `solve!` not called yet"""
        Convex.MAXWIDTH[] = old_maxwidth

        # solved problem
        x = Variable()
        p = satisfy(x >= 0)
        output = solve!(p, SCS.Optimizer(verbose=0, eps=1e-6))
        @test sprint(show, p) == """
                minimize
                └─ 0
                subject to
                └─ >= constraint (affine)
                   ├─ real variable ($(Convex.show_id(x)))
                   └─ 0

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
        @test length(Convex.show_id(x)) == length("id: ") + length(string(x.id_hash))
        Convex.MAXDIGITS[] = old_maxdigits

    end

    @testset "vartype and set_vartype" begin
        for x in (Variable(), Variable(1), ComplexVariable(2, 2))
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
    end

    @testset "Constructors" begin

        # Constructors with sign
        for sgn in (Positive(), NoSign())
            for x in    [   # tuple size
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
                            Variable(sgn, :Bin),  ]
                @test x isa Variable
                @test x isa Convex.AbstractVariable
                @test sign(x) == sgn
                @test x.sign == sgn
            end
        end

        # constructors without sign
        for x in    [   # tuple size
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
                        Variable(:Bin),  ]
            @test x isa Variable
            @test x isa Convex.AbstractVariable
            @test sign(x) == NoSign()
            @test x.sign == NoSign()

            Convex.sign!(x, Positive())
            @test sign(x) == Positive()
            @test x.sign == Positive()
        end

        # ComplexVariable
        for x in    [   # tuple size
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
                    Variable(ComplexSign()),  ]
            @test x isa Variable
            @test x isa Convex.AbstractVariable
            @test sign(x) == ComplexSign()
            @test x.sign == ComplexSign()
        end
        
        for vt in (BinVar, IntVar), V in (ComplexVariable, Semidefinite, HermitianSemidefinite)
            @test_throws Any V(2; vartype=vt)
        end

        for vt in (:Bin, :Int), V in (Semidefinite, HermitianSemidefinite, ComplexVariable)
            @test_throws Any V(2, vt)
        end

        # Semidefinite
        for x in [
                Variable((2,2), :Semidefinite),
                Variable(2,2, :Semidefinite),
                ComplexVariable((2,2), :Semidefinite),
                ComplexVariable(2,2, :Semidefinite),
                HermitianSemidefinite((2,2)),
                HermitianSemidefinite(2, 2),
                HermitianSemidefinite(2),
                Semidefinite((2,2)),
                Semidefinite(2, 2),
                Semidefinite(2),
            ]
            @test length(constraints(x)) == 1
            @test constraints(x)[] isa Convex.SDPConstraint
        end

        @test_throws ErrorException HermitianSemidefinite(2,3)
        @test_throws ErrorException Semidefinite(2,3)
    end

    @testset "ConicObj" for T = [UInt32, UInt64]
        c = ConicObj()
        z = zero(T)
        @test !haskey(c, z)
        c[z] = (1, 1)
        @test c[z] == (1, 1)
        x = T[]
        for (k, v) in c
            push!(x, k)
        end
        @test x == collect(keys(c))
        d = copy(c)
        @test d !== c
    end

    @testset "length and size" begin
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
    end

    @testset "lastindex and axes" begin
        x = Variable(2, 3)
        @test axes(x) == (Base.OneTo(2), Base.OneTo(3))
        @test axes(x, 1) == Base.OneTo(2)
        @test lastindex(x) == 6
        @test lastindex(x, 2) == 3
        y = x[:,end]
        @test y isa AbstractExpr
        @test size(y) == (2, 1)
    end

    @testset "Parametric constants" begin
        z = Constant([1.0 0.0im; 0.0 1.0])
        @test z isa Constant{Matrix{Complex{Float64}}}

        # Helper functions
        @test Convex.ispos(1)
        @test Convex.ispos(0)
        @test !Convex.ispos(-1)
        @test Convex.ispos([0,1,0])
        @test !Convex.ispos([0,-1,0])
        @test Convex.isneg(-1)
        @test Convex.isneg(0)
        @test !Convex.isneg(1)
        @test Convex.isneg([0,-1,0])
        @test !Convex.isneg([0,1,0])
        @test Convex._size(3) == (1, 1)
        @test Convex._sign(3) == Positive()
        @test Convex._size([-1,1,1]) == (3, 1)
        @test Convex._sign([-1,1,1]) == NoSign()
        @test Convex._sign([-1,-1,-1]) == Negative()
        @test Convex._size([0 0; 0 0]) == (2, 2)
        @test Convex._sign([0 0; 0 0]) == Positive()
        @test Convex._size(0 + 1im) == (1, 1)
        @test Convex._sign(0 + 1im) == ComplexSign()

        @test Convex.imag_conic_form(Constant(1.0)) == [0.0]
        @test Convex.imag_conic_form(Constant([1.0, 2.0])) == [0.0, 0.0]
    end

    @testset "#341: Evaluate for constants" begin
        A = rand(4,4)
        @test evaluate(Constant(A)) ≈ copy(A)
        @test Constant(A).size == (4,4)
        b = rand(4)
        @test evaluate(Constant(b)) ≈ copy(b)
        @test Constant(b).size == (4,1)
        c = 1.0
        @test evaluate(Constant(c)) ≈ c
        @test Constant(c).size == (1,1)

        @test evaluate(sumlargesteigs(Variable(4, 4), 0)) == 0
        @test evaluate(sumlargest(Variable(4), 0)) == 0
        @test evaluate(sumsmallest(Variable(4), 0)) == 0
    end

    @testset "Base.vect" begin
    # Issue #223: ensure we can make vectors of variables
        @test size([Variable(2), Variable(3, 4)]) == (2,)
    end

    @testset "Iteration" begin
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
    end
    # returns [21]; not sure why
    # context("iteration") do
    #     x = Variable(2,3)
    #     s = sum([xi for xi in x])
    #     x.value = [1 2 3; 4 5 6]
    #     @fact evaluate(s) --> 21
    # end

    @testset "DCP warnings" begin
        # default is to log
        @test_logs (:warn, r"not DCP compliant") Convex.NotDcp()
        
        @eval Convex.emit_dcp_warnings() = false
        @test_logs  Convex.NotDcp()
        @eval Convex.emit_dcp_warnings() = true
        @test_logs (:warn, r"not DCP compliant") Convex.NotDcp()

    end

    @testset "`add_constraints!` (#380)" begin
        x = Variable(3, 3)
        p = minimize(norm_1(x))
        y = randn(3, 3)
        c = (norm2(x-y) < 1)
        @test length(p.constraints) == 0
        add_constraint!(p, c)
        @test length(p.constraints) == 1
        empty!(p.constraints)
        add_constraints!(p, c)
        @test length(p.constraints) == 1
        empty!(p.constraints)
        c2 = (norm2(x-rand(3,3)) < 3)
        add_constraints!(p, [c, c2])
        @test length(p.constraints) == 2
    end

    @testset "`diagm` (#401)" begin
        x = Variable(3)
        @test diagm(x) isa AbstractExpr
    end
end
