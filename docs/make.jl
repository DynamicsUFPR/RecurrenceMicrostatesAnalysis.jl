##
##      Library Documenter.jl to build the documentation page.
using Documenter, RMA

##
##      Build docs.
makedocs(
    sitename = "RMA.jl",
    format = Documenter.HTML(
        prettyurls = true
    ),
    pages = [
        "Welcome" => "index.md",
        "Guide" => [
            "Quick Start" => "quickstart.md",
            "Theoretical Overview" => "theory.md",
            "Distributions" => "distributions.md",
            "RQA" => "rqa.md",
            "Utils" => "utils.md",
            "Performance Tips" => "performance.md",
        ],
        "Reference" => [
            "Public API" => "api.md",
            "Motifs: shapes and sampling" => "motifs.md",
            "Recurrence functions" => "recurrence.md"
        ],
        "Examples" => "examples.md",
        "Bibliography" => "bib.md",
        "Developers and Researchers" => "dev.md",
        "Release notes" => "release_notes.md",
    ]
)