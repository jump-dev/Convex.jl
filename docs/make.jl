using Documenter, Convex, Literate, Pkg

# Needed to run GR headless on Travis
previous_GKSwstype = get(ENV, "GKSwstype", "")
ENV["GKSwstype"] = "100"

@info "Building examples..."

filename(str) = first(splitext(last(splitdir(str))))
filename_to_name(str) = uppercasefirst(replace(replace(filename(str), "-" => " "), "_" => " "))

fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")

literate_path = joinpath(@__DIR__(), "examples_literate")
build_path =  joinpath(@__DIR__, "src", "examples")
rm(build_path; force=true, recursive=true)
isdir(build_path) || mkdir(build_path)


@info "[Examples] Preparing notebooks..."

notebooks_path = joinpath(@__DIR__, "notebooks")
rm(notebooks_path, recursive=true, force=true)
mkdir(notebooks_path)

for dir in readdir(literate_path)
    dir_path = joinpath(literate_path, dir)
    isdir(dir_path) || continue
    @info "Processing directory $dir"
    notebook_dir = joinpath(notebooks_path, dir)
    isdir(notebook_dir) || mkdir(notebook_dir)
    for file in readdir(dir_path)
        file_path = joinpath(dir_path, file)
            out_path = joinpath(notebooks_path, dir, file)
        if endswith(file, ".jl")
            Literate.notebook(file_path, notebook_dir, execute=false)
        else
            cp(file_path, out_path)
        end
    end
end

# Copy `Project.toml` to notebooks
cp(joinpath(@__DIR__, "Project.toml"), joinpath(notebooks_path, "Project.toml"))

# Add a README file to notebooks
open(joinpath(notebooks_path, "README.md"), "w") do io
    print(io, """
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
    """)
end

# zip up the notebooks directory
zip_path = joinpath(build_path, "notebooks.zip")
run(Cmd(`zip $zip_path -r notebooks`; dir = @__DIR__))

@info "[Examples] Preparing markdown files..."

for dir in readdir(literate_path)
    dir_path = joinpath(literate_path, dir)
    isdir(dir_path) || continue
    @info "Processing directory $dir"
    build_dir = joinpath(build_path, dir)
    isdir(build_dir) || mkdir(build_dir)
    for file in readdir(dir_path)
        file_path = joinpath(dir_path, file)
            out_path = joinpath(build_path, dir, file)
        if endswith(file, ".jl")
            postprocess = function(content)
                """
                All of the examples can be found in Jupyter notebook form [here](../$(filename(zip_path)).zip).
                """ * content
            end
            Literate.markdown(file_path, build_dir; preprocess = fix_math_md, documenter = true, postprocess =  postprocess)
        else
            cp(file_path, out_path)
        end
    end
end

@info "Starting `makedocs`"

# Build nav tree for examples
function nav_dir(dir, path)
    sort([ joinpath("examples", dir, file) for file in readdir(path) if endswith(file, ".md") && file != "index.md" ])
end

examples_nav = [ filename_to_name(dir) => nav_dir(dir, joinpath(build_path, dir)) for dir in readdir(build_path) if isdir(joinpath(build_path, dir)) ]

makedocs(;
    modules = [Convex],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Quick Tutorial" => "quick_tutorial.md",
        "Basic Types" => "types.md",
        "Supported Operations" => "operations.md",
        "Complex-domain Optimization" => "complex-domain_optimization.md",
        "Solvers" => "solvers.md",
        "FAQ" => "faq.md",
        "Advanced" => "advanced.md",
        "Problem Depot" => "problem_depot.md",
        "Contributing" => "contributing.md",
        "Credits" => "credits.md",
        "Examples" => examples_nav,
    ],
    repo = "https://github.com/JuliaOpt/Convex.jl/blob/{commit}{path}#L{line}",
    sitename = "Convex.jl")

deploydocs(repo = "github.com/JuliaOpt/Convex.jl.git")

# restore the environmental variable `GKSwstype`.
ENV["GKSwstype"] = previous_GKSwstype;
