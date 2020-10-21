seed(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.(seed(flux,v,cfg)...)

seed_position(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.(seed_position(flux,v,cfg)...)

seed_strength(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.(seed_strength(flux,v,cfg)...)

seed(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.(seed(flux,v,cfg)...,Elements.blobradius(v))

seed_position(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.(seed_position(flux,v,cfg)...,Elements.blobradius(v))

seed_strength(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.(seed_strength(flux,v,cfg)...,Elements.blobradius(v))
