
seed(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.(seed(circulation,v,cfg)...)

seed_position(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.(seed_position(circulation,v,cfg)...)

seed_strength(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.(seed_strength(circulation,v,cfg)...)

seed(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.(seed(circulation,v,cfg)...,Elements.blobradius(v))

seed_position(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.(seed_position(circulation,v,cfg)...,Elements.blobradius(v))

seed_strength(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.(seed_strength(circulation,v,cfg)...,Elements.blobradius(v))
