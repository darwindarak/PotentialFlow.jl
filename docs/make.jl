using Documenter, VortexModel

makedocs(
    format =:html,
    sitename = "VortexModel.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [#"manual/quickstart.md",
                     "manual/elements.md",
                     "manual/velocities.md",
                     "manual/timemarching.md",
                     "manual/noflowthrough.md"
                     ],
        "Internals" => [ "internals/properties.md"]
    ],
    assets = ["assets/custom.css"],
    strict = true
)

if "DOCUMENTER_KEY" in keys(ENV)
    deploydocs(
     repo = "github.com/darwindarak/VortexModel.jl.git",
     target = "build",
     deps = nothing,
     make = nothing,
     julia = "0.6"
    )
end
