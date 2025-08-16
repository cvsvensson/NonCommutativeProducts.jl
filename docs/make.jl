using NonCommutativeProducts
using Documenter

DocMeta.setdocmeta!(NonCommutativeProducts, :DocTestSetup, :(using NonCommutativeProducts); recursive=true)

makedocs(;
    modules=[NonCommutativeProducts],
    authors="Viktor Svensson",
    sitename="NonCommutativeProducts.jl",
    format=Documenter.HTML(;
        canonical="https://cvsvensson.github.io/NonCommutativeProducts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cvsvensson/NonCommutativeProducts.jl",
    devbranch="main",
)
