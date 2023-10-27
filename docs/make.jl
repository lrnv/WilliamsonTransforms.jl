using WilliamsonTransforms
using Documenter

DocMeta.setdocmeta!(WilliamsonTransforms, :DocTestSetup, :(using WilliamsonTransforms); recursive=true)

makedocs(;
    modules=[WilliamsonTransforms],
    authors="Oskar Laverny <oskar.laverny@gmail.com> and contributors",
    repo="https://github.com/lrnv/WilliamsonTransforms.jl/blob/{commit}{path}#{line}",
    sitename="WilliamsonTransforms.jl",
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
    repo="github.com/lrnv/WilliamsonTransforms.jl",
    devbranch="master",
)
