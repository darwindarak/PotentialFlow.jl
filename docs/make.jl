using Documenter, PotentialFlow

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

makedocs(
    sitename = "PotentialFlow.jl",
    doctest = true,
    clean = true,
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
    #format = Documenter.HTML(assets = ["assets/custom.css"]),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
#   strict = true
)

#if "DOCUMENTER_KEY" in keys(ENV)
    #deploydocs(;
    # repo = "github.com/darwindarak/PotentialFlow.jl.git",
    #)
deploydocs(
    repo = "github.com/darwindarak/PotentialFlow.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
    #versions = "v^"
    )
#end
