using Documenter, Express

makedocs(;
    modules=[Express],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/MineralsCloud/Express.jl/blob/{commit}{path}#L{line}",
    sitename="Express.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/MineralsCloud/Express.jl",
)
