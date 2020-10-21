for f in (:seed, :seed_position, :seed_strength)
  @eval $f(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.($f(flux,v,cfg)...)
end

for f in (:seed, :seed_position, :seed_strength)
  @eval $f(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.($f(flux,v,cfg)...,Elements.blobradius(v))
end
