using NonCommutativeProducts
using Documenter
using DocumenterCodeBlocks

DocMeta.setdocmeta!(NonCommutativeProducts, :DocTestSetup, :(using NonCommutativeProducts); recursive=true)

makedocs(;
    modules=[NonCommutativeProducts],
    authors="Viktor Svensson, William Samuelson",
    sitename="NonCommutativeProducts.jl",
    repo=Documenter.Remotes.GitHub("cvsvensson", "NonCommutativeProducts.jl"),
    format=Documenter.HTML(;
        canonical="https://cvsvensson.github.io/NonCommutativeProducts.jl",
        edit_link="main",
    ),
    plugins=[CodeBlocks()],
    checkdocs=:public,
    pages=[
        "Home" => "index.md",
        "Multiple types" => "multiple_types.md",
        "Traversing expressions" => "traversal.md",
        "Linear algebra" => "linear_algebra.md",
        "Gradings" => "gradings.md",
        "Performance tips" => "performance.md",
        "API reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/cvsvensson/NonCommutativeProducts.jl",
    devbranch="main",
)
