using Documenter, PotentialFlow

makedocs(
    format =:html,
    sitename = "PotentialFlow.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/quickstart.md",
                     "manual/elements.md",
                     "manual/velocities.md",
                     "manual/timemarching.md",
                     "manual/noflowthrough.md",
                     "manual/motions.md"
                     ],
        "Internals" => [ "internals/properties.md"]
    ],
    assets = ["assets/custom.css"],
    strict = true,
    doctest = !is_windows()
)


if "DOCUMENTER_KEY" in keys(ENV)
    deploydocs(
     repo = "github.com/darwindarak/PotentialFlow.jl.git",
     target = "build",
     deps = nothing,
     make = nothing,
     julia = "0.6"
    )
end
