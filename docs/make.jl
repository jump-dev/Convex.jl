using Documenter, Convex

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
