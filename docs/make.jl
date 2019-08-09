using Documenter, Convex

makedocs(;
            modules = [Convex],
            format=Documenter.HTML(),
            # pages = [
                # "Home" => "index.md",
            # ],
            repo="https://github.com/JuliaOpt/Convex.jl/blob/{commit}{path}#L{line}",
            sitename="Convex.jl",
            authors="")
