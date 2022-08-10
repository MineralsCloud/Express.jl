using ExpressBase
using Documenter

DocMeta.setdocmeta!(ExpressBase, :DocTestSetup, :(using ExpressBase); recursive=true)

makedocs(;
    modules=[ExpressBase],
    authors="Reno <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/ExpressBase.jl/blob/{commit}{path}#{line}",
    sitename="ExpressBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation guide" => "installation.md",
            "Contributing" => "contributing.md",
        ],
        "Library" => "public.md",
        "Troubleshooting" => "troubleshooting.md",
    ],
)
