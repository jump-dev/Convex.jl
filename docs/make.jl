using Documenter, Convex, Literate, Pkg

# Needed to run GR headless on Travis
previous_GKSwstype = get(ENV, "GKSwstype", "")
ENV["GKSwstype"] = "100"

@info "Building examples..."

filename(str) = first(splitext(last(splitdir(str))))
filename_to_name(str) = uppercasefirst(replace(replace(filename(str), "-" => " "), "_" => " "))

function fix_math_md(content)
    replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
end

literate_path = joinpath(@__DIR__(), "..", "examples", "literate")
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
                ```@setup $(filename(file_path))
                PREVIOUS_PROJECT = Base.active_project();
                using Pkg
                Pkg.activate(@__DIR__)
                Pkg.instantiate()
                ```
                All of the examples can be found in Jupyter notebook form [here](../$(filename(zip_path)).zip).
                """ * content * """
                ```@setup $(filename(file_path))
                Pkg.activate(PREVIOUS_PROJECT)
                ```
                """
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

# Get current environment, since our examples will change environments
current_project = Base.active_project()

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

# Restore environment
Pkg.activate(current_project)

deploydocs(repo = "github.com/JuliaOpt/Convex.jl.git")


ENV["GKSwstype"] = previous_GKSwstype
