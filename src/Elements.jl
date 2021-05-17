module Elements

export Element, Singleton, Group, kind, @kind, circulation, flux, streamfunction,
       complexpotential, conftransform, inverse_conftransform, jacobian, image,
       blobradius

using ..Properties

import ..Utils: ComplexDual, Dual, ComplexGradientConfig, seed!, seed,
                vector_mode_gradient, vector_mode_jacobian, checktag,chunksize,
                chunk_mode_gradient_expr, Partials, valtype, extract_gradient_chunk!

abstract type Element end

#== Trait definitions ==#

abstract type Singleton end
abstract type Group end

function kind end

macro kind(element, k)
    esc(quote
        Elements.kind(::$element) = $k
        Elements.kind(::Type{$element}) = $k
        end)
end

@kind ComplexF64 Singleton

kind(::ComplexDual) = Singleton
kind(::Dual) = Singleton
kind(::AbstractArray{T}) where {T <: Union{Element, ComplexF64, Dual, ComplexDual}} = Group
kind(::Tuple) = Group


# Convenience functions to define wrapper types
# e.g. vortex sheets as a wrapper around vortex blobs
unwrap_src(e) = unwrap(e)
unwrap_targ(e) = unwrap(e)
unwrap(e) = e

## Actual Definitions of Properties

"""
    Elements.position(src::Element)

Returns the complex position of a potential flow element.
This is a required method for all `Element` types.

# Example

```jldoctest
julia> point = Vortex.Point(1.0 + 0.0im, 1.0);

julia> Elements.position(point)
1.0 + 0.0im

julia> points = Vortex.Point.([1.0im, 2.0im], 1.0);

julia> Elements.position.(points)
2-element Array{Complex{Float64},1}:
 0.0 + 1.0im
 0.0 + 2.0im
```
"""
@property begin
    signature = position(src::Source)
    stype = ComplexF64
end

"""
    Elements.blobradius(src)

Returns the blob radius of a regularized potential flow element.
If the element is singular, returns zero.
"""
@property begin
    signature = blobradius(src::Source)
    stype = Float64
end

"""
    Elements.circulation(src)

Returns the total circulation contained in `src`.

# Example

```jldoctest
julia> points = Vortex.Point.([1.0im, 2.0im], [1.0, 2.0]);

julia> blobs = Vortex.Blob.([1.0im, 2.0im], [1.0, 2.0], 0.1);

julia> Elements.circulation(points[1])
1.0

julia> Elements.circulation(points)
3.0

julia> Elements.circulation((points, blobs))
6.0

julia> Elements.circulation.(points)
2-element Array{Float64,1}:
 1.0
 2.0

julia> Elements.circulation.((points, blobs))
(3.0, 3.0)

julia> Elements.circulation(Source.Point(rand(), rand()))
0.0

julia> Elements.circulation(Source.Blob(rand(), rand(), rand()))
0.0
```
"""
@property begin
    signature = circulation(src::Source)
    reduce = (+)
    stype = Float64
end

"""
    Elements.flux(src)

Returns the flux through a unit circle induced by `src`.

# Example

```jldoctest
julia> points = Source.Point.([1.0im, 2.0im], [1.0, 2.0]);

julia> blobs = Source.Blob.([1.0im, 2.0im], [1.0, 2.0], 0.1);

julia> Elements.flux(points[1])
1.0

julia> Elements.flux((points, blobs))
6.0

julia> Elements.flux.(points)
2-element Array{Float64,1}:
 1.0
 2.0

julia> Elements.flux.((points, blobs))
(3.0, 3.0)

julia> Elements.flux(Vortex.Point(rand(), rand()))
0.0

julia> Elements.flux(Vortex.Blob(rand(), rand(), rand()))
0.0
```
"""
@property begin
    signature = flux(src::Source)
    reduce = (+)
    stype = Float64
end

@doc raw"""
    Elements.impulse(src)

Return the aerodynamic impulse of `src` about (0,0):
```math
P := \int \boldsymbol{x} \times \boldsymbol{\omega}\,\mathrm{d}A.
```
This is a required method for all vortex types.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(1.0im, -π, 0.1));

julia> Elements.impulse(sys[1])
3.141592653589793 + 0.0im

julia> Elements.impulse(sys)
0.0 + 0.0im
```
"""
@property begin
    signature = impulse(src::Source)
    reduce = (+)
    stype = ComplexF64
end

@doc raw"""
    Elements.angularimpulse(src)

Return the aerodynamic angular impulse of `src` about (0,0):
```math
\Pi := \frac{1}{2}\int \boldsymbol{x} \times (\boldsymbol{x}\times\boldsymbol{\omega})\,\mathrm{d}A.
```
This is a required method for all vortex types.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(1.0im, -π, 0.1));

julia> Elements.angularimpulse(sys[1])
3.141592653589793 + 0.0im

julia> Elements.angularimpulse(sys)
0.0 + 0.0im
```
"""
@property begin
    signature = angularimpulse(src::Source)
    reduce = (+)
    stype = ComplexF64
end

@property begin
    signature = streamfunction(targ::Target, src::Source)
    preallocator = allocate_streamfunction
    stype = Float64
end

@property begin
    signature = complexpotential(targ::Target, src::Source)
    preallocator = allocate_complexpotential
    stype = ComplexF64
end

@doc raw"""
    Elements.conftransform(src,body)

Return the conformally transformed position of `src` via the transform defined
by `body`. In this context, the position(s) in `src` is interpreted as
coordinates in the circle plane.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(2.0im, -π, 0.1));

julia> b = ConformalBody([1/4,0,1/4]);

julia> Elements.conftransform(sys,b)
(0.0 + 0.0im, 0.0 + 0.375im)
```
"""
@property begin
    signature = conftransform(src::Source,b)
    preallocator = allocate_conftransform
    stype = ComplexF64
end

@doc raw"""
    Elements.inverse_conftransform(src,body)

Return the inverse conformally transformed position of `src` via the transform defined
by `body`. In this context, the position(s) in `src` is interpreted as
coordinates in the physical plane.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(2.0im, -π, 0.1));

julia> b = ConformalBody([1/4,0,1/4]);

julia> Elements.inverse_conftransform(sys,b)
(0.0 + 0.0im, 0.0 + 0.375im)
```
"""
@property begin
    signature = inverse_conftransform(src::Source,b)
    preallocator = allocate_inv_conftransform
    stype = ComplexF64
end


@doc raw"""
    Elements.jacobian(src,body)

Return the Jacobian of the conformal transform at the position of `src` via the transform defined
by `body`. In this context, the position(s) in `src` is interpreted as
coordinates in the circle plane.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(2.0im, -π, 0.1));

julia> b = ConformalBody([1/4,0,1/4]);

julia> Elements.jacobian(sys,b)
(0.5 + 0.0im, 0.0 + 0.5im)

```
"""
@property begin
    signature = jacobian(src::Source,b)
    preallocator = allocate_jacobian
    stype = ComplexF64
end


@doc raw"""
    Elements.image(src,body)

Return the image position of `src` in the circle-plane representation of `body`.
In this context, the position(s) in `src` is interpreted as
coordinates in the circle plane.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(2.0im, -π, 0.1));

julia> b = PowerBody([1/4,0,1/4],zero(ComplexF64),0.0)

julia> Elements.image(sys,b)
(0.0 + 1.0im, 0.0 + 0.5im)

```
"""
@property begin
    signature = image(src::Source,b)
    preallocator = allocate_image
    stype = ComplexF64
end




"""
    seed_strength(T,v::Vector{Element},i::Int)

Given a collection `v` of points or blobs, create a copy of the collection with
the `i`th element's strength replaced by a unit dual. The entire output
collection has positions of `Dual` type, but only the replaced element has unit
partials; the others have partials equal to zero. `T` is a `Tag`, and can be set
to `Nothing`.
"""
function seed_strength end

"""
    property_type(::Element)

Return the type of the properties of the element.
"""
function property_type end

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
  seed!(strduals,complex(strfunc.(v)),index)
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
  seed!(posduals,position(v),index)
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




end
