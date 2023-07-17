using ExpressBase
using Documenter

DocMeta.setdocmeta!(
    ExpressBase, :DocTestSetup, :(using ExpressBase, ExpressBase.Files); recursive=true
)

makedocs(;
    modules=[ExpressBase],
    authors="singularitti <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/ExpressBase.jl/blob/{commit}{path}#{line}",
    sitename="ExpressBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/ExpressBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => ["Installation guide" => "installation.md"],
        "API Reference" => "public.md",
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/ExpressBase.jl",
    devbranch="main",
)
