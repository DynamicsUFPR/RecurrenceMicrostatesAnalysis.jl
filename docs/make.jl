##
##      Library Documenter.jl to build the documentation page.
using Documenter
using DocumenterCitations

include("../src/RecurrenceMicrostatesAnalysis.jl")
using .RecurrenceMicrostatesAnalysis

bib = CitationBibliography("reference.bib")

##
##      Build docs.
makedocs(
    sitename = "RecurrenceMicrostatesAnalysis.jl",
    format = Documenter.HTML(
        prettyurls = true
    ),
    modules = [RecurrenceMicrostatesAnalysis],
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
        # "Developers and Researchers" => "dev.md",
        # "Release notes" => "release_notes.md",
    ],
    plugins = [bib]
)

deploydocs(
    repo = "github.com/DynamicsUFPR/RecurrenceMicrostatesAnalysis.jl.git",
    branch = "docs",
    versions = [
        "dev" => "dev",
        "stable" => "main",
        "v0.2.24" => "v0.2.24",
    ]
)