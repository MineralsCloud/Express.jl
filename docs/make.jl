using Express
using Documenter

makedocs(;
    modules=[Express],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/Express.jl/blob/{commit}{path}#L{line}",
    sitename="Express.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/Express.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => ["Installation" => "install.md", "Development" => "develop.md"],
        "Examples" => ["hcp-GaN example" => "examples/GaN.md"],
        "API by module" => [],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/Express.jl",
)
