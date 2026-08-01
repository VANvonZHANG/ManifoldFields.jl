using Documenter
using ManifoldFields

makedocs(
    sitename = "ManifoldFields.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    modules = [ManifoldFields],
    pages = [
        "Home" => "index.md",
        "Roadmap" => "roadmap.md",
        "API Reference" => "api.md"
    ],
    warnonly = [:missing_docs]
)

if get(ENV, "MANIFOLDFIELDS_DOCS_DEPLOY", "false") == "true"
    deploydocs(
        repo = "github.com/VANvonZHANG/ManifoldFields.jl.git",
        devbranch = "main",
        push_preview = true
    )
end
