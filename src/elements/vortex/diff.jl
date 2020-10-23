for f in (:seed_zeros, :seed_position, :seed_strength)
  @eval $f(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.($f(circulation,v,cfg)...)
end

for f in (:seed_zeros, :seed_position, :seed_strength)
  @eval $f(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.($f(circulation,v,cfg)...,Elements.blobradius(v))
end
