__precompile__()

module PotentialFlow

using Reexport

include("Properties.jl")
include("Utils.jl")

include("Elements.jl")
include("Motions.jl")

@reexport using .Elements
@reexport using .Motions

include("TimeMarching.jl")
@reexport using .TimeMarching

#== Some Built-in Potential Flow Elements ==#

export Vortex, Source, Sheets, Plates, Plate, Bodies, PowerBody

include("elements/Points.jl")
include("elements/Blobs.jl")
include("elements/Sheets.jl")

include("elements/Vortex.jl")
include("elements/Source.jl")

include("elements/Plates.jl")
include("elements/Bodies.jl")

import .Plates: Plate
import .Bodies: PowerBody


#== Plot Recipes ==#

include("plot_recipes.jl")


end
