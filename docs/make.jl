using Convex
import Documenter
import Literate

const SKIP_EXAMPLES = get(ENV, "CONVEX_SKIP_EXAMPLES", false) == "true"
const BUILD_PATH = joinpath(@__DIR__, "src", "examples")
const LITERATE_PATH = joinpath(@__DIR__(), "examples_literate")
const NOTEBOOKS_PATH = joinpath(@__DIR__, "notebooks")

rm(BUILD_PATH; force = true, recursive = true)
if !isdir(BUILD_PATH)
    mkdir(BUILD_PATH)
end

filename(str) = first(splitext(last(splitdir(str))))

function filename_to_name(str)
    return uppercasefirst(
        replace(replace(filename(str), "-" => " "), "_" => " "),
    )
end

fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")

examples_nav = Any[]

if SKIP_EXAMPLES
    @info "Skipping examples"
else
    @info "Building examples..."
    @info "[Examples] Preparing notebooks..."
    rm(NOTEBOOKS_PATH, recursive = true, force = true)
    mkdir(NOTEBOOKS_PATH)
    for dir in readdir(LITERATE_PATH)
        dir_path = joinpath(LITERATE_PATH, dir)
        if !isdir(dir_path)
            continue
        end
        @info "Processing directory $dir"
        notebook_dir = joinpath(NOTEBOOKS_PATH, dir)
        if !isdir(notebook_dir)
            mkdir(notebook_dir)
        end
        for file in readdir(dir_path)
            file_path = joinpath(dir_path, file)
            out_path = joinpath(NOTEBOOKS_PATH, dir, file)
            if endswith(file, ".jl")
                Literate.notebook(file_path, notebook_dir, execute = false)
            else
                cp(file_path, out_path)
            end
        end
    end
    # Copy `Project.toml` to notebooks
    cp(
        joinpath(@__DIR__, "Project.toml"),
        joinpath(NOTEBOOKS_PATH, "Project.toml"),
    )
    # Add a README file to notebooks
    open(joinpath(NOTEBOOKS_PATH, "README.md"), "w") do io
        print(
            io,
            """
  # Convex.jl example notebooks

  Start Julia in this directory and set the project flag to point to this directory. E.g. run the command
  ```julia
  julia --project=.
  ```
  in this directory.

  Then add `IJulia` if it's not installed already in your global environment by
  ```julia
  pkg> add IJulia
  ```

  Also call `instantiate` to download the required packages:
  ```julia
  pkg> instantiate
  ```

  Then launch Jupyter:
  ```julia
  julia> using IJulia

  julia> notebook(dir=pwd(); detached=true)
  ```

  This should allow you to try any of the notebooks.
  """,
        )
        return
    end
    # zip up the notebooks directory
    zip_path = joinpath(BUILD_PATH, "notebooks.zip")
    run(Cmd(`zip $zip_path -r notebooks`; dir = @__DIR__))
    @info "[Examples] Preparing markdown files..."
    for dir in readdir(LITERATE_PATH)
        dir_path = joinpath(LITERATE_PATH, dir)
        if !isdir(dir_path)
            continue
        end
        @info "Processing directory $dir"
        build_dir = joinpath(BUILD_PATH, dir)
        if !isdir(build_dir)
            mkdir(build_dir)
        end
        for file in readdir(dir_path)
            file_path = joinpath(dir_path, file)
            out_path = joinpath(BUILD_PATH, dir, file)
            if endswith(file, ".jl")
                postprocess =
                    (content) -> begin
                        block_name = replace(filename(file), r"\s+" => "_")
                        return """
                               All of the examples can be found in Jupyter notebook form [here](../$(filename(zip_path)).zip).

                               ```@setup $(block_name)
                               __START_TIME = time_ns()
                               @info "Starting example $(filename(file))"
                               ```
                               """ *
                               content *
                               """
               ```@setup $(block_name)
               __END_TIME = time_ns()
               elapsed = string(round((__END_TIME - __START_TIME)*1e-9; sigdigits = 3), "s")
               @info "Finished example $(filename(file)) after " * elapsed
               ```
               """
                    end
                Literate.markdown(
                    file_path,
                    build_dir;
                    preprocess = fix_math_md,
                    documenter = true,
                    postprocess = postprocess,
                )
            else
                cp(file_path, out_path)
            end
        end
    end
    # Build nav tree for examples
    function nav_dir(dir, path)
        return sort([
            joinpath("examples", dir, file) for
            file in readdir(path) if endswith(file, ".md") && file != "index.md"
        ])
    end
    for dir in readdir(BUILD_PATH)
        path = joinpath(BUILD_PATH, dir)
        if isdir(path)
            push!(examples_nav, filename_to_name(dir) => nav_dir(dir, path))
        end
    end
end

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
#  Build annd release
# ==============================================================================

Documenter.makedocs(
    sitename = "Convex.jl",
    repo = "https://github.com/jump-dev/Convex.jl/blob/{commit}{path}#L{line}",
    # TODO(odow): uncomment this once all docstrings are in the manual
    # modules = [Convex],
    format = Documenter.HTML(;
        ansicolor = true,
        collapselevel = 1,
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    strict = true,
    pages = [
        "Introduction" => [
            "Home" => "index.md",
            "introduction/installation.md",
            "introduction/quick_tutorial.md",
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
