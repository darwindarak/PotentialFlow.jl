module VortexModel
using Reexport

include("Vortex.jl")

@reexport using .Vortex

export Vortex

end
