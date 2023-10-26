using WilliamsonTransform
using Documenter

DocMeta.setdocmeta!(WilliamsonTransform, :DocTestSetup, :(using WilliamsonTransform); recursive=true)

makedocs(;
    modules=[WilliamsonTransform],
    authors="Oskar Laverny <oskar.laverny@gmail.com> and contributors",
    repo="https://github.com/lrnv/WilliamsonTransform.jl/blob/{commit}{path}#{line}",
    sitename="WilliamsonTransform.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lrnv.github.io/WilliamsonsTranform.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lrnv/WilliamsonTransform.jl",
    devbranch="master",
)
