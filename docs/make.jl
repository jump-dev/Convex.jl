using Documenter, Convex, Literate, Pkg

previous_GKSwstype = get(ENV, "GKSwstype", "")
ENV["GKSwstype"] = "100"

@info "Building examples"

filename(str) = first(splitext(last(splitdir(str))))
filename_to_name(str) = uppercasefirst(replace(replace(filename(str), "-" => " "), "_" => " "))

function fix_math_md(content)
    replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
end


literate_path = joinpath(@__DIR__(), "..", "examples", "literate")
build_path =  joinpath(@__DIR__, "src", "examples")
rm(build_path; force=true, recursive=true)
isdir(build_path) || mkdir(build_path)
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
                """ * content * """
                ```@setup $(filename(file_path))
                Pkg.activate(PREVIOUS_PROJECT)
                ```
                """
            end
            documenter = true
            Literate.markdown(file_path, build_dir; preprocess = fix_math_md, documenter = documenter, postprocess =  documenter ? postprocess : identity)
            Literate.notebook(file_path, build_dir; execute=false)
            name = string(filename_to_name(file), " (webpage)")
            path = filename(file) * ".md"
        else
            cp(file_path, out_path)
        end
    end
end


function nav_dir(dir, path)
    sort([ joinpath("examples", dir, file) for file in readdir(path) if endswith(file, ".md") && file != "index.md" ])
end

examples_nav = [ filename_to_name(dir) => nav_dir(dir, joinpath(build_path, dir)) for dir in readdir(build_path) if isdir(joinpath(build_path, dir)) ]

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

Pkg.activate(@__DIR__)

deploydocs(repo = "github.com/JuliaOpt/Convex.jl.git")

ENV["GKSwstype"] = previous_GKSwstype
