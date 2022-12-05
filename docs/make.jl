using RedClust, Documenter, Literate

# Parse the basic example file
function preprocess_md(content)
    content = replace(content, "#setup!" => "# ```@setup @__NAME__\n#=") # hacky workaround to create @setup blocks in Literate
    content = replace(content, "#!setup" => "=#\n# ```")
    return content
end
function preprocess_jl(content)
    content = replace(content, "#setup!" => "")
    content = replace(content, "#!setup" => "")
    content = replace(content, r" *nothing *# *hide *(\r\n|\n)" => "") # remove the nothing # hide lines when writing to .jl files
    return content
end
inputdir = joinpath(@__DIR__, "..", "examples")
inputfile = joinpath(inputdir, "basic_example.jl")
outputdir = joinpath(@__DIR__, "src", "_generated")
# Create the example file in the docs
Literate.markdown(inputfile, outputdir; 
name = "basic_example", preprocess = preprocess_md)
# Create the actual example file
Literate.script(inputfile, outputdir; 
name = "basic_example", preprocess = preprocess_jl, keep_comments=true)

# Build the documentation HTML pages
Documenter.Writers.HTMLWriter.HTML(ansicolor = true)
makedocs(;
    modules=[RedClust],
    authors="Abhinav Natarajan <abhinav.v.natarajan@gmail.com>",
    repo="https://github.com/abhinavnatarajan/RedClust.jl/blob/{commit}{path}#{line}",
    sitename="RedClust.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://abhinavnatarajan.github.io/RedClust.jl",
        edit_link="main",
        assets=String[],
        ansicolor=true
    ),
    pages=[
    "Introduction" => "index.md",
    "Example" => joinpath("_generated", "basic_example.md"),
    "Reference" => "reference.md",
    "Changelog" => "changelog.md"
    ]
)

deploydocs(;
    repo="github.com/abhinavnatarajan/RedClust.jl",
    devbranch="master"
)