ComplexGradientConfig(f::F,v::Vector{<:Element},w...) where {F} =
      ComplexGradientConfig(f,position(v),w...)


"""
    seed_zeros(v::Vector{Element},cfg::ComplexGradientConfig)

Given a collection `v` of points or blobs, create a copy of the collection with
the element's positions and strengths replaced by complex duals with zeros.
`cfg` specifies the type and size of duals. The use of complex duals (in spite of the real-valued
strength) ensures correct treatment for gradient calculation, since some
calculations will be complex-valued.
"""
function seed_zeros(strength::F,v::Vector{<:Element},cfg::ComplexGradientConfig) where F
  posduals = convert.(eltype(cfg.duals),position(v))
  strduals = convert.(eltype(cfg.duals),complex(strength.(v)))
  return posduals, real(strduals)
end

function seed!(posduals::AbstractArray{<:ComplexDual{T,V,M}},strduals::AbstractArray{<:ComplexDual{T,V,M}},strfunc::F,v::Vector{<:Element}) where {F,T,V,M}
  seed!(posduals,position(v))
  seed!(strduals,complex(strfunc.(v)))
  return posduals, real(strduals)
end

function seed!(posduals::AbstractArray{<:ComplexDual{T,V,M}},strduals::AbstractArray{<:ComplexDual{T,V,M}},strfunc::F,v::Vector{<:Element},index) where {F,T,V,M}
  seed!(posduals,position(v),index)
  seed!(strduals,complex(strfunc.(v)),index)
  return posduals, real(strduals)
end

"""
    seed_position(v::Vector{Element},cfg::ComplexGradientConfig)

Given a collection `v` of points or blobs, create a copy of the collection with
the element's positions replaced by complex duals with unit seeds and with
the element's strengths replaced by complex duals with zeros. `cfg` specifies
the type and size of duals. The length of the partials in these duals is twice
the number of elements, to enable computation of gradients with respect to
all positions.
"""
function seed_position!(posduals::AbstractArray{<:ComplexDual{T,V,M}},strduals::AbstractArray{<:ComplexDual{T,V,M}},
                        strfunc::F,v::Vector{<:Element},rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}}) where {T,F,N,M,V}
  seed!(posduals,position(v),rseeds,iseeds)
  seed!(strduals,complex(strfunc.(v)))
  return posduals, real(strduals)
end

function seed_position!(posduals::AbstractArray{<:ComplexDual{T,V,M}},strduals::AbstractArray{<:ComplexDual{T,V,M}},
                        strfunc::F,v::Vector{<:Element},index,rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}},
                        chunksize=N) where {T,F,N,M,V}
  seed!(posduals,position(v),index,rseeds,iseeds,chunksize)
  seed!(strduals,complex(strfunc.(v)),index,zero(Partials{M,V}),zero(Partials{M,V}),chunksize)
  return posduals, real(strduals)
end

seed_position(strfunc::F,v::Vector{<:Element},cfg::ComplexGradientConfig) where F =
    seed_position!(copy(cfg.duals),copy(cfg.duals),strfunc,v,cfg.rseeds,cfg.iseeds)


"""
    seed_strength(v::Vector{Element},cfg::ComplexGradientConfig)

Given a collection `v` of points or blobs, create a copy of the collection with
the element's strengths replaced by complex duals with unit seeds and with
the element's positions replaced by complex duals with zeros. `cfg` specifies
the type and size of duals. The use of complex duals (in spite of the real-valued
strength) ensures correct treatment for gradient calculation, since some
calculations will be complex-valued.
"""
function seed_strength!(posduals::AbstractArray{<:ComplexDual{T,V,M}},strduals::AbstractArray{<:ComplexDual{T,V,M}},
                        strfunc::F,v::Vector{<:Element},rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}}) where {T,F,N,M,V}
  seed!(posduals,position(v))
  seed!(strduals,complex(strfunc.(v)),rseeds,iseeds)
  return posduals, real(strduals)
end

function seed_strength!(posduals::AbstractArray{<:ComplexDual{T,V,M}},strduals::AbstractArray{<:ComplexDual{T,V,M}},
                        strfunc::F,v::Vector{<:Element},index,
                        rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}},chunksize=N) where {T,F,N,M,V}
  seed!(posduals,position(v),index,zero(Partials{M,V}),zero(Partials{M,V}),chunksize)
  seed!(strduals,complex(strfunc.(v)),index,rseeds,iseeds,chunksize)
  return posduals, real(strduals)
end



seed_strength(strfunc::F,v::Vector{<:Element},cfg::ComplexGradientConfig) where F =
    seed_strength!(copy(cfg.duals),copy(cfg.duals),strfunc,v,cfg.rseeds,cfg.iseeds)


## Gradient

function gradient_position(f,v::Vector{<:Element},cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, v)) where {T}
    checktag(T, f, property_type(eltype(v)))
    if chunksize(cfg) == length(v)
      vector_mode_gradient(f,v,cfg,vector_mode_dual_eval_position)
    else
      chunk_mode_gradient_position(f,v,cfg)
    end
end

function gradient_strength(f,v::Vector{<:Element},cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, v)) where {T}
    checktag(T, f, property_type(eltype(v)))
    if chunksize(cfg) == length(v)
      dz, dzstar = vector_mode_gradient(f,v,cfg,vector_mode_dual_eval_strength)
    else
      dz, dzstar = chunk_mode_gradient_strength(f,v,cfg)
    end
    return dz .+ dzstar
end

@inline vector_mode_dual_eval_position(f::F, v::Vector{<:Element},
          cfg::ComplexGradientConfig) where {F} = f(seed_position(v,cfg))

@inline vector_mode_dual_eval_strength(f::F, v::Vector{<:Element},
          cfg::ComplexGradientConfig) where {F} = f(seed_strength(v,cfg))

@inline vector_mode_dual_eval_param(f::F, v::Tuple{Vector{<:Element},T},
          cfg::ComplexGradientConfig) where {F,T} = f(seed_zeros(v[1],cfg),seed(v[2],cfg))

function jacobian_position(f,v::Vector{<:Element},cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, v)) where {T}
    checktag(T, f, property_type(eltype(v)))
    if chunksize(cfg) == length(v)
      vector_mode_jacobian(f,v,cfg,vector_mode_dual_eval_position)
    else
      chunk_mode_jacobian_position(f,v,cfg)
    end
end

function jacobian_strength(f,v::Vector{<:Element},cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, v)) where {T}
    checktag(T, f, property_type(eltype(v)))
    if chunksize(cfg) == length(v)
      dz, dzstar = vector_mode_jacobian(f,v,cfg,vector_mode_dual_eval_strength)
    else
      dz, dzstar = chunk_mode_jacobian_strength(f,v,cfg)
    end
    return dz .+ dzstar
end

function jacobian_param(f,v::Tuple{Vector{<:Element},S},
        cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, ComplexF64[v[2]])) where {S,T}
    checktag(T, f, v[2])
    dz, dzstar = vector_mode_jacobian(f,v,cfg,vector_mode_dual_eval_param)
    return real(dz .+ dzstar)
end

for ftype = (:position, :strength)
  fname = Symbol("chunk_mode_gradient_",ftype)
  sname = Symbol("seed_",ftype,"!")
  @eval function $fname(f::F, z::Vector{<:Element}, cfg::ComplexGradientConfig{T,V,N}) where {F,T,V,N}
     $(chunk_mode_gradient_expr(quote
                                  posduals = copy(cfg.duals)
                                  strduals = copy(cfg.duals)
                                  zdual = seed_zeros(z,cfg)
                                  seed!(zdual,z,posduals,strduals)
                                end,
                                quote
                                  dz = similar(z, Complex{valtype(ydual)})
                                  dzstar = similar(z, Complex{valtype(ydual)})
                                end,
                                :(dz, dzstar),
                                :(ydual = f(zdual)),
                                :(),
                                :($sname(zdual, z, posduals, strduals, index, rseeds, iseeds, chunksize)),
                                :(seed!(zdual, z, posduals, strduals, index)),
                                :(extract_gradient_chunk!(T, dz, dzstar, ydual, index, chunksize)),
                                :()))
  end
end

for ftype = (:position, :strength)
  fname = Symbol("chunk_mode_jacobian_",ftype)
  sname = Symbol("seed_",ftype,"!")
  @eval function $fname(f::F, z::Vector{<:Element}, cfg::ComplexGradientConfig{T,V,N}) where {F,T,V,N}
     $(chunk_mode_gradient_expr(quote
                                  posduals = copy(cfg.duals)
                                  strduals = copy(cfg.duals)
                                  zdual = seed_zeros(z,cfg)
                                  seed!(zdual,z,posduals,strduals)
                                end,
                                quote
                                  dz = similar(z, Complex{valtype(ydual)},length(ydual), zlen)
                                  dzstar = similar(z, Complex{valtype(ydual)},length(ydual), zlen)
                                  dz_reshaped = reshape_jacobian(dz, ydual, zdual)
                                  dzstar_reshaped = reshape_jacobian(dzstar, ydual, zdual)
                                end,
                                :(dz, dzstar),
                                :(ydual = f(zdual)),
                                :(ydual isa AbstractArray || throw(JACOBIAN_ERROR)),
                                :($sname(zdual, z, posduals, strduals, index, rseeds, iseeds, chunksize)),
                                :(seed!(zdual, z, posduals, strduals, index)),
                                :(extract_jacobian_chunk!(T, dz_reshaped, dzstar_reshaped, ydual, index, chunksize)),
                                :()))
  end
end
