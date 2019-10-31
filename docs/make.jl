using Documenter, Convex

@info "Building example notebooks"
# build examples notebooks
notebooks_path = joinpath(@__DIR__, "src", "notebooks")
include(joinpath(@__DIR__, "..", "examples", "build.jl"))

@info "Generating `examples.md`"
# generate the `examples.md` page with links to the notebooks
filename_to_name(str) = uppercasefirst(replace(replace(splitext(str)[1], "-" => " "), "_" => " "))
open(joinpath(@__DIR__, "src", "examples.md"), "w") do io
    println(io, "# Examples")
    for example_dir in readdir(notebooks_path)
        name = filename_to_name(example_dir)
        println(io, "\n## $name\n")
        for notebook in readdir(joinpath(notebooks_path, example_dir))
            endswith(notebook, ".ipynb") || continue
            path = joinpath("notebooks", example_dir, notebook)
            name = string(filename_to_name(notebook), " (Jupyter notebook)")
            println(io, "* [$name]($path)")
        end
    end
end

@info "Starting `makedocs`"

makedocs(;
    modules = [Convex],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Quick Tutorial" => "quick_tutorial.md",
        "Basic Types" => "types.md",
        "Supported Operations" => "operations.md",
        "Examples" => "examples.md",
        "Complex-domain Optimization" => "complex-domain_optimization.md",
        "Solvers" => "solvers.md",
        "FAQ" => "faq.md",
        "Advanced" => "advanced.md",
        "Problem Depot" => "problem_depot.md",
        "Contributing" => "contributing.md",
        "Credits" => "credits.md"
    ],
    repo = "https://github.com/JuliaOpt/Convex.jl/blob/{commit}{path}#L{line}",
    sitename = "Convex.jl")


deploydocs(repo = "github.com/JuliaOpt/Convex.jl.git")
