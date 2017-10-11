module Elements

export Element, Singleton, Group, kind, @kind, circulation

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

@kind Complex128 Singleton
kind(::AbstractArray{T}) where {T <: Union{Element, Complex128}} = Group
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
    stype = Complex128
end

"""
    Vortex.circulation(src)

Returns the total circulation contained in `src`.

# Example

```jldoctest
julia> points = Vortex.Point.([1.0im, 2.0im], [1.0, 2.0]);

julia> Elements.circulation(points[1])
1.0

julia> Elements.circulation(points)
3.0

julia> Elements.circulation.(points)
2-element Array{Float64,1}:
 1.0
 2.0
```
"""
@property begin
    signature = circulation(src::Source)
    reduce = (+)
    stype = Float64
end

doc"""
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
    stype = Complex128
end

end
