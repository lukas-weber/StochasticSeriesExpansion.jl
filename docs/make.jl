using Documenter, StochasticSeriesExpansion

makedocs(
    sitename = "StochasticSeriesExpansion.jl",
    format = Documenter.HTML(prettyurls = false),
    checkdocs = :all,
    pages = [
        "index.md",
        # "tutorial.md",
        "Implementing custom models" => ["interfaces.md", "sse_data.md"],
    ],
)

deploydocs(repo = "github.com/lukas-weber/StochasticSeriesExpansion.jl.git")
