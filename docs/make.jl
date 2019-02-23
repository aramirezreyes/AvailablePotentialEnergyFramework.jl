using Documenter, SAMtools

makedocs(
    modules = [SAMtools],
    format = :html,
    checkdocs = :exports,
    sitename = "SAMtools.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/aramirezreyes/SAMtools.jl.git",
)
