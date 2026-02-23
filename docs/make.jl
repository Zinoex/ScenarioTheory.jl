using ScenarioTheory
using Documenter

DocMeta.setdocmeta!(ScenarioTheory, :DocTestSetup, :(using ScenarioTheory); recursive=true)

makedocs(;
    modules=[ScenarioTheory],
    authors="Frederik Baymler Mathiesen <frederik@baymler.com> and contributors",
    sitename="ScenarioTheory.jl",
    format=Documenter.HTML(;
        canonical="https://Zinoex.github.io/ScenarioTheory.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Zinoex/ScenarioTheory.jl",
    devbranch="main",
)
