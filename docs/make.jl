push!(LOAD_PATH,"../src/")

using Documenter
using BoxCoxTransformations

makedocs(
    sitename = "BoxCoxTransformations.jl",
    doctest = false,
    authors = "Oliver Lindemann",
    pages = [
        "index.md",
    ],
)

deploydocs(;
    repo = "github.com/lindemann09/BoxCoxTransformations.jl",
    push_preview = true,
    devbranch="main"
)