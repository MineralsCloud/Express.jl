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
            "Installation" => "install.md",
            "Development" => "develop.md",
            "Running" => [
                "Configuration files" => "configuration.md",
                "How to run `Express` from command line" => "run.md",
            ],
        ],
        "Troubleshooting" => "troubleshooting.md",
        "Other questions" => "questions.md",
        "API by module" => ["`EquationOfStateWorkflow` module" => "api/EquationOfStateWorkflow.md"],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/Express.jl",
)
