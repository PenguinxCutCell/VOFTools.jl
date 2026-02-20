using Documenter
using VOFTools

makedocs(
    modules = [VOFTools],
    authors = "Fastaxx and contributors",
    repo = "https://github.com/PenguinxCutCell/VOFTools.jl/blob/{commit}{path}#{line}",
    sitename = "VOFTools.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/VOFTools.jl",
        repolink = "https://github.com/PenguinxCutCell/VOFTools.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Geometry Types" => "geometry-types.md",
        "2D Tools" => "tools2d.md",
        "3D Tools" => "tools3d.md",
        "Mesh Generators" => "mesh-generators.md",
        "Reference" => "95-reference.md",
    ],
    pagesonly = true,
    warnonly = true,
)

deploydocs(
    repo = "github.com/PenguinxCutCell/VOFTools.jl",
    push_preview = true,
)
