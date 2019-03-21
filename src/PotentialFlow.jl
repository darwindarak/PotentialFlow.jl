__precompile__()

module PotentialFlow

using Reexport

include("Properties.jl")
include("Utils.jl")

include("Elements.jl")
include("RigidBodyMotions.jl")
include("Motions.jl")


@reexport using .Elements
@reexport using .Motions
@reexport using .RigidBodyMotions

include("TimeMarching.jl")
@reexport using .TimeMarching

#== Some Built-in Potential Flow Elements ==#

export Vortex, Source, Sheets, Plates, Plate, Freestreams, Freestream,
       Doublets, Doublet, Bodies

include("elements/Points.jl")
include("elements/Blobs.jl")
include("elements/Sheets.jl")

include("elements/Vortex.jl")
include("elements/Source.jl")

include("elements/Doublets.jl")
include("elements/Freestream.jl")

include("elements/Plates.jl")
include("elements/Bodies.jl")

@reexport using SchwarzChristoffel

import .Plates: Plate
import .Bodies: ConformalBody





#== Plot Recipes ==#

include("plot_recipes.jl")


end
