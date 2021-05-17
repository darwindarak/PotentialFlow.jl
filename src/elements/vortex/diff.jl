for f in (:seed_zeros, :seed_position, :seed_strength)
  @eval $f(v::Vector{<:Point},cfg::ComplexGradientConfig) =
      Point.($f(circulation,v,cfg)...)

  @eval $f(v::Vector{<:Blob},cfg::ComplexGradientConfig) =
      Blob.($f(circulation,v,cfg)...,Elements.blobradius(v))

end

for fiip in (:seed_position!, :seed_strength!)
  @eval function $fiip(vseed::Vector{<:Point},v::Vector{<:Point},posduals,strduals,
                      index,rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}},
                      chunksize=N) where {N,V,M}
      $fiip(posduals,strduals,circulation,v,index,rseeds,iseeds,chunksize)
      vseed .= Point.(posduals,real(strduals))
  end

  @eval function $fiip(vseed::Vector{<:Blob},v::Vector{<:Blob},posduals,strduals,
                      index,rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}},
                      chunksize=N) where {N,V,M}
      $fiip(posduals,strduals,circulation,v,index,rseeds,iseeds,chunksize)
      vseed .= Blob.(posduals,real(strduals),Elements.blobradius(v))
  end
end

for fiip in (:seed!,)
  @eval function $fiip(vseed::Vector{<:Point},v::Vector{<:Point},posduals,strduals)
      $fiip(posduals,strduals,circulation,v)
      vseed .= Point.(posduals,real(strduals))
  end

  @eval function $fiip(vseed::Vector{<:Blob},v::Vector{<:Blob},posduals,strduals)
      $fiip(posduals,strduals,circulation,v)
      vseed .= Blob.(posduals,real(strduals),Elements.blobradius(v))
  end

  @eval function $fiip(vseed::Vector{<:Point},v::Vector{<:Point},posduals,strduals,index) where {N,V,M}
      $fiip(posduals,strduals,circulation,v,index)
      vseed .= Point.(posduals,real(strduals))
  end

  @eval function $fiip(vseed::Vector{<:Blob},v::Vector{<:Blob},posduals,strduals,index) where {N,V,M}
      $fiip(posduals,strduals,circulation,v,index)
      vseed .= Blob.(posduals,real(strduals),Elements.blobradius(v))
  end
end
