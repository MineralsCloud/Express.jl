using Express
using Documenter

DocMeta.setdocmeta!(Express, :DocTestSetup, :(using Express); recursive=true)

makedocs(;
    modules=[Express],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/Express.jl/blob/{commit}{path}#{line}",
    sitename="Express.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/Express.jl",
        assets=String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Installation guide" => "installation.md",
            "Contributing" => "contributing.md",
            "Running" => [
                "Configuration files" => "run/configuration.md",
                "How to run `Express` from command line" => "run/cli.md",
                "Tracking and monitoring jobs in a workflow" => "run/jobs.md",
            ],
        ],
        "Troubleshooting" => "troubleshooting.md",
        "FAQ" => "faq.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/Express.jl",
)
