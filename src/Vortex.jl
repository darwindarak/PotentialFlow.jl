"""
The `Vortex` module provides basic vortex elements:

- [`Vortex.Point`](@ref)
- [`Vortex.Blob`](@ref)
- [`Vortex.Sheet`](@ref)
- [`Vortex.Plate`](@ref)

and functions that act on collections vortex element types.

## Exported functions

- [`induce_velocity`](@ref)
- [`induce_velocity!`](@ref)
- [`mutually_induce_velocity!`](@ref)
- [`self_induce_velocity!`](@ref)
- [`advect!`](@ref)
- [`allocate_velocity`](@ref)
- [`reset_velocity!`](@ref)
- [`@get`](@ref)

## Useful unexported functions

- [`Vortex.position`](@ref)
- [`Vortex.circulation`](@ref)
- [`Vortex.impulse`](@ref)
- [`Vortex.unit_impulse`](@ref)

"""
module Vortex

export allocate_velocity, reset_velocity!,
    induce_velocity, induce_velocity!,
    mutually_induce_velocity!, self_induce_velocity!,
    advect!, @get

include("Utils.jl")
using .Utils

abstract type Element end
abstract type PointSource <: Element end
abstract type CompositeSource <: Element end

const Collection = Union{AbstractArray, Tuple}
const PointArray = AbstractArray{T} where {T <: PointSource}
const TargetTypes = (Complex128, PointSource, CompositeSource, Collection)

"""
    Vortex.position(src::PointSource)

Returns the complex position of a `PointSource` type vortex element
This is a required method for all subtypes of `PointSource`.
"""
function position end

"""
    Vortex.circulation(src)

Returns the total circulation contained in `src`
This is a required method for all vortex types.
"""
circulation(vs::Collection) = mapreduce(circulation, +, 0.0, vs)

"""
    Vortex.impulse(src)

Returns the aerodynamic impulse of `src` relative to (0, 0)
This is a required method for all vortex types.
"""
impulse(vs::Collection) = mapreduce(impulse, +, 0.0, vs)

"""
    allocate_velocity(srcs)

Allocate arrays of `Complex128` to match the structure of `srcs`

For example:

```julia
julia> points = Vortex.Point.(rand(Complex128, 2), rand(2));
julia> blobs  = Vortex.Blob.(rand(Complex128, 3), rand(3), rand(3));
julia> allocate_velocity(points)
(Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im])
```
"""
function allocate_velocity(array::AbstractArray{T}) where {T <: Union{PointSource, Complex128}}
    zeros(Complex128, length(array))
end

allocate_velocity(group::Tuple) = map(allocate_velocity, group)

"""
    reset_velocity!(vels[, srcs])

Set all velocities in `vels` to zero

If `srcs` is provided, then the arrays in `vels` are resized their source counterpart, if necessary.
"""
reset_velocity!(array::Vector{Complex128}) = fill!(array, zero(Complex128))
reset_velocity!(::Void, ::Any) = nothing

function reset_velocity!(group::Tuple)
    for ẋ in group
        reset_velocity!(ẋ)
    end
    group
end

function reset_velocity!(array::Vector{Complex128}, v)
    if length(array) != length(v)
        resize!(array, length(v))
    end
    fill!(array, zero(Complex128))
end

function reset_velocity!(group::T, src) where {T <: Tuple}
    for i in 1:length(group)
        reset_velocity!(group[i], src[i])
    end
    group
end

"""
    induce_velocity(target, element)

Compute the velocity induced by `element` on `target`

`target` can be:

- a `Complex128`
- a subtype of `Vortex.PointSource`
- an array or tuple of vortex elements

while the `element` can be:

- any subtype of `Vortex.Element`
- an array or tuple of vortex elements
"""
function induce_velocity(target::PointSource, source)
    induce_velocity(Vortex.position(target), source)
end

function induce_velocity(z::Complex128, sources::Collection)
    w = zero(Complex128)
    for source in sources
        w += induce_velocity(z, source)
    end
    w
end

function induce_velocity(targets::Collection, source)
    ws = allocate_velocity(targets)
    induce_velocity!(ws, targets, source)
end

"""
    induce_velocity!(vels, target, element)

Compute the velocity induced by `element` on `target` and store the result in `vels`

`vels` should be the output of a call to [`allocate_velocity`](@ref),
`target` can be an array or tuple of vortex elements,
while the `element` can be:
- any subtype of `Vortex.Element`
- an array or tuple of vortex elements
"""
function induce_velocity!(ws::AbstractArray, targets::AbstractArray, source)
    for i in 1:length(targets)
        ws[i] += induce_velocity(targets[i], source)
    end
    ws
end

induce_velocity!(::Void, ::Any, ::Any) = nothing

function induce_velocity!(ws::Tuple, targets::Tuple, source)
    for i in 1:length(targets)
        induce_velocity!(ws[i], targets[i], source)
    end
    ws
end

"""
    mutually_induce_velocity!(vs₁, vs₂, e₁, e₂)

Compute the mutually induced velocities between `e₁` and `e₂` and
store the results in `vs₁` and `vs₂`

The default implementation simply calls [`induce_velocity!`](@ref) twice.
This method is meant to be overwritten to take advantage of symmetries in certain pairwise vortex interations.
For example, the velocity kernel for a point vortex is antisymmetric,
so in computing the mutually induced velocities of two arrays of point vortices,
we can half the number of calls to the velocity kernel.
"""
function mutually_induce_velocity!(ws₁, ws₂, v₁, v₂)
    induce_velocity!(ws₁, v₁, v₂)
    induce_velocity!(ws₂, v₂, v₁)
    nothing
end

"""
    self_induce_velocity!(vels, elements)

Compute the self induced velocity of one or more vortex elements

This involves a recursive call to `self_induce_velocity!` and pairwise calls to [`mutually_induce_velocity!`](@ref).
"""
function self_induce_velocity!(ws, group::Tuple)
    N = length(group)
    for s in 1:N
        self_induce_velocity!(ws[s], group[s])
        for t in s+1:N
            mutually_induce_velocity!(ws[s], ws[t], group[s], group[t])
        end
    end
    nothing
end

"""
    self_induce_velocity(src)

Compute the self induced velocity of a point source.
If this method is not implemented for a subtype of `Vortex.PointSource`, the fallback function returns `0.0+0.0im`
"""
self_induce_velocity(::PointSource) = zero(Complex128)

function self_induce_velocity!(ws, array::PointArray)
    N = length(array)
    for s in 1:N
        ws[s] += self_induce_velocity(array[s])
        for t in s+1:N
            ws[s] += induce_velocity(array[s], array[t])
            ws[t] += induce_velocity(array[t], array[s])
        end
    end
    nothing
end

"""
    advect(src::PointSource, velocity, Δt)

Return a new vortex element that represents `src` advected by `velocity` over `Δt`
If this method is implemented by any type `T <: PointSource`,
then an array of type `AbstractArray{T}` can be passed in the first two arguments of [`advect!`](@ref).
"""
function advect end

function advect!{T <: Collection}(vs₊::T, vs₋::T, w, Δt)
    for (i, v) in enumerate(vs₋)
        vs₊[i] = advect(v, w[i], Δt)
    end
    nothing
end

"""
    advect!(srcs₊, srcs₋, vels, Δt)

Moves the elements in `srcs₋` by their corresponding velocity in `vels` over the interval `Δt` and store the results in `src₊`
`srcs₋` and `srcs₊` can be either a array of vortex elements or a tuple.
"""
function advect!(group₊::T, group₋::T, ws, Δt) where {T <: Tuple}
    for i in 1:length(group₋)
        advect!(group₊[i], group₋[i], ws[i], Δt)
    end
    nothing
end

@submodule "elements/Points"
@submodule "elements/Blobs"
@submodule "elements/Sheets"
@submodule "elements/Plates"

end
