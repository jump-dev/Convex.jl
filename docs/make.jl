# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

using Convex
import Documenter
import Literate
import Test

const SKIP_EXAMPLES = get(ENV, "CONVEX_SKIP_EXAMPLES", false) == "true"

function _file_list(full_dir, relative_dir, extension)
    return map(
        file -> joinpath(relative_dir, file),
        filter(file -> endswith(file, extension), sort(readdir(full_dir))),
    )
end

function _literate_directory(dir)
    for filename in _file_list(dir, dir, ".md")
        rm(filename)
    end
    for filename in _file_list(dir, dir, ".jl")
        if endswith(filename, "antidiag.jl")
            continue
        elseif endswith(filename, "data.jl")
            continue
        end
        # `include` the file to test it before `#src` lines are removed. It is
        # in a testset to isolate local variables between files.
        Test.@testset "$(filename)" begin
            # we can't use just `mod = Module()`
            # as some examples use `include`.
            # TODO- remove `include`'s from examples
            mod = @eval module $(gensym()) end
            Base.include(mod, filename)
        end
        Literate.markdown(filename, dir; documenter = true)
    end
    return
end

if !SKIP_EXAMPLES
    Test.@testset "Examples" verbose = true begin
        for (root, dir, files) in walkdir(joinpath(@__DIR__, "src", "examples"))
            _literate_directory.(joinpath.(root, dir))
        end
    end
end

root = joinpath(@__DIR__, "src", "examples")
examples_nav = Any[
    uppercasefirst(replace(dir, "_" => " ")) => [
        joinpath("examples", dir, file) for
        file in readdir(joinpath(root, dir)) if endswith(file, ".md")
    ] for dir in readdir(root) if isdir(joinpath(root, dir))
]

# ==============================================================================
#  Modify the release notes
# ==============================================================================

function fix_release_line(
    line::String,
    url::String = "https://github.com/jump-dev/Convex.jl",
)
    # (#XXXX) -> ([#XXXX](url/issue/XXXX))
    while (m = match(r"\(\#([0-9]+)\)", line)) !== nothing
        id = m.captures[1]
        line = replace(line, m.match => "([#$id]($url/issues/$id))")
    end
    # ## vX.Y.Z -> [vX.Y.Z](url/releases/tag/vX.Y.Z)
    while (m = match(r"\#\# (v[0-9]+.[0-9]+.[0-9]+)", line)) !== nothing
        tag = m.captures[1]
        line = replace(line, m.match => "## [$tag]($url/releases/tag/$tag)")
    end
    return line
end

open(joinpath(@__DIR__, "src", "changelog.md"), "r") do in_io
    open(joinpath(@__DIR__, "src", "release_notes.md"), "w") do out_io
        for line in readlines(in_io; keep = true)
            write(out_io, fix_release_line(line))
        end
    end
end

# ==============================================================================
#  Build and release
# ==============================================================================

Documenter.makedocs(
    sitename = "Convex.jl",
    # TODO(odow): uncomment this once all docstrings are in the manual
    # modules = [Convex],
    format = Documenter.HTML(;
        ansicolor = true,
        collapselevel = 1,
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    warnonly = [:cross_references],
    pages = [
        "Introduction" => [
            "Home" => "index.md",
            "introduction/installation.md",
            "introduction/quick_tutorial.md",
            "introduction/dcp.md",
            "introduction/faq.md",
        ],
        "Examples" => examples_nav,
        "Manual" => [
            "manual/types.md",
            "manual/operations.md",
            "Complex-domain Optimization" => "manual/complex-domain_optimization.md",
            "manual/solvers.md",
            "manual/advanced.md",
        ],
        "Developer Docs" => [
            "developer/problem_depot.md",
            "developer/contributing.md",
            "developer/credits.md",
        ],
        "reference.md",
        "release_notes.md",
    ],
)

Documenter.deploydocs(
    repo = "github.com/jump-dev/Convex.jl.git";
    push_preview = true,
)
