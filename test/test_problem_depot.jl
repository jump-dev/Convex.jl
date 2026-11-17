# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

module TestProblemDepot

using Convex
using Test

import Hypatia
# import ECOS
import GLPK
import MathOptInterface as MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function _test_problem_depot_load(T)
    Convex.ProblemDepot.foreach_problem() do name, func
        @testset "$name" begin
            func(Val(false), 0.0, 0.0, T) do problem
                @test problem isa Convex.Problem{T} # check numeric type
                model =
                    () -> MOI.Utilities.MockOptimizer(MOI.Utilities.Model{T}())
                context = Convex.Context(problem, model)
                @test sprint(show, context) isa AbstractString
                @test sprint(show, problem) isa AbstractString
                return
            end
        end
    end
    return
end

test_problem_depot_load_Float64() = _test_problem_depot_load(Float64)
test_problem_depot_load_BigFloat() = _test_problem_depot_load(BigFloat)

function test_Hypatia_warmstarts()
    Convex.ProblemDepot.run_tests(; exclude = [r"mip"]) do p
        ret = solve!(p, Hypatia.Optimizer; silent = true, warmstart = true)
        # verify post-solve printing does not error
        @test sprint(show, p) isa AbstractString
        return ret
    end
    return
end

function test_Hypatia()
    Convex.ProblemDepot.run_tests(; exclude = [r"mip"]) do p
        return solve!(p, Hypatia.Optimizer; silent = true)
    end
    return
end

# Disabling ECOS due to this segfault:
# https://github.com/jump-dev/ECOS.jl/issues/144
# function test_ecos()
#     # For `rational_norm` problems:
#     # >   MathOptInterface.UnsupportedConstraint{MathOptInterface.VectorAffineFunction{Float64}, MathOptInterface.NormCone}: `MathOptInterface.VectorAffineFunction{Float64}`-in-`MathOptInterface.NormCone` constraint is not supported by the model.
#     # Missing a bridge?
#     Convex.ProblemDepot.run_tests(; exclude = [r"mip", r"sdp", r"rational_norm"]) do p
#         return solve!(p, ECOS.Optimizer; silent = true)
#     end
#     return
# end

function test_GLPK()
    # Note this is an inclusion, not exclusion;
    # we only test GLPK with MIPs, since those we can't test elsewhere.
    Convex.ProblemDepot.run_tests([r"mip"]) do p
        return solve!(p, GLPK.Optimizer; silent = true)
    end
    return
end

end

TestProblemDepot.runtests()
