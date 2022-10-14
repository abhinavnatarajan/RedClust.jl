using RedClust
using Documenter

DocMeta.setdocmeta!(RedClust, :DocTestSetup, :(using RedClust); recursive=true)

makedocs(;
    modules=[RedClust],
    authors="Abhinav Natarajan <abhinav.v.natarajan@gmail.com>",
    repo="https://github.com/ab      hinavnatarajan/RedClust.jl/blob/{commit}{path}#{line}",
    sitename="RedClust.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://abhinavnatarajan.github.io/RedClust.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
    "Introduction" => "index.md",
    "Reference" => "reference.md",
    "Changelog" => "changelog.md"
    ]
)

deploydocs(;
    repo="github.com/abhinavnatarajan/RedClust.jl",
    devbranch="master"
)