module Elements

export Element, Singleton, Group, kind, @kind, circulation, flux, streamfunction,
       conftransform, inverse_conftransform, jacobian, image

using ..Properties

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
kind(::AbstractArray{T}) where {T <: Union{Element, ComplexF64}} = Group
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
 0.0+1.0im
 0.0+2.0im
```
"""
@property begin
    signature = position(src::Source)
    stype = ComplexF64
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

raw"""
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

raw"""
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

raw"""
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

raw"""
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


raw"""
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


raw"""
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

end
