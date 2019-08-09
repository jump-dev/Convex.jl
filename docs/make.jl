using Documenter, Convex

makedocs(;
            modules = [Convex],
            format=Documenter.HTML(),
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
                "Optimizing in a Loop" => "loop.md",
                "Advanced" => "advanced.md",
                "Contributing" => "contributing.md",
                "Credits" => "credits.md"
            ],
            repo="https://github.com/JuliaOpt/Convex.jl/blob/{commit}{path}#L{line}",
            sitename="Convex.jl")


deploydocs(;
    repo="https://github.com/JuliaOpt/Convex.jl",
)

