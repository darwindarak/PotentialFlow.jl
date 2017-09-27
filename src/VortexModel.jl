__precompile__()

module VortexModel
using Reexport

include("Vortex.jl")
include("TimeMarching.jl")

@reexport using .Vortex
@reexport using .TimeMarching

export Vortex

end
