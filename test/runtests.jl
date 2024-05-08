# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

using Test, Aqua, Convex

# Seed random number stream to improve test reliability
import Random
Random.seed!(2)

@testset "Convex" begin
    include("test_utilities.jl")
    include("test_atoms.jl")
    include("test_constraints.jl")
    include("test_problem_depot.jl")
    include("MOI_wrapper.jl")
end

@testset "Aqua" begin
    Aqua.test_all(Convex; ambiguities=false, piracies=false)
    # Convex currently has some light piracy of `hcat`, `vcat` and `hvcat`
    Aqua.test_piracies(Convex; treat_as_own=[hcat, vcat, hvcat])
    # Convex currently has some ambiguities
    Aqua.test_ambiguities(Convex; broken=true)
end
