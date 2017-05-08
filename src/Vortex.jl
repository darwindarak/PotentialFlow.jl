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
    advect, advect!, @get

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

# Example

```jldoctest
julia> point = Vortex.Point(1.0 + 0.0im, 1.0);

julia> Vortex.position(point)
1.0 + 0.0im

julia> points = Vortex.Point.([1.0im, 2.0im], 1.0);

julia> Vortex.position.(points)
2-element Array{Complex{Float64},1}:
 0.0+1.0im
 0.0+2.0im
```
"""
function position end

"""
    Vortex.circulation(src)

Returns the total circulation contained in `src`
This is a required method for all vortex types.

# Example

```jldoctest
julia> points = Vortex.Point.([1.0im, 2.0im], [1.0, 2.0]);

julia> Vortex.circulation(points[1])
1.0

julia> Vortex.circulation(points)
3.0

julia> Vortex.circulation.(points)
2-element Array{Float64,1}:
 1.0
 2.0
```
"""
circulation(vs::Collection) = mapreduce(circulation, +, 0.0, vs)

doc"""
    Vortex.impulse(src)

Return the aerodynamic impulse of `src` about (0,0):
```math
P := \int \boldsymbol{x} \times \boldsymbol{\omega}\,\mathrm{d}A.
```
This is a required method for all vortex types.

# Example

```jldoctest
julia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(1.0im, -π, 0.1));

julia> Vortex.impulse(sys[1])
3.141592653589793 + 0.0im

julia> Vortex.impulse(sys)
0.0 + 0.0im
```
"""
impulse(vs::Collection) = mapreduce(impulse, +, 0.0, vs)

"""
    allocate_velocity(srcs)

Allocate arrays of `Complex128` to match the structure of `srcs`

# Example

```jldoctest
julia> points = Vortex.Point.(rand(Complex128, 2), rand(2));

julia> blobs  = Vortex.Blob.(rand(Complex128, 3), rand(3), rand(3));

julia> allocate_velocity(points)
2-element Array{Complex{Float64},1}:
 0.0+0.0im
 0.0+0.0im

julia> allocate_velocity((points, blobs))
(Complex{Float64}[0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im])
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

# Example

```jldoctest
julia> ẋs = (rand(Complex128, 1), rand(Complex128, 1))
(Complex{Float64}[0.236033+0.346517im], Complex{Float64}[0.312707+0.00790928im])

julia> points = Vortex.Point.(rand(Complex128, 2), rand(2));

julia> blobs  = Vortex.Blob.(rand(Complex128, 3), rand(3), rand(3));

julia> reset_velocity!(ẋs, (points, blobs));

julia> ẋs
(Complex{Float64}[0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im])
```
"""
reset_velocity!(array::Vector{Complex128}) = fill!(array, zero(Complex128))

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

# Example
```jldoctest
julia> z = rand(Complex128)
0.23603334566204692 + 0.34651701419196046im

julia> point = Vortex.Point(z, rand());

julia> srcs = Vortex.Point.(rand(Complex128, 10), rand(10));

julia> induce_velocity(z, srcs[1])
0.08722212007570912 + 0.14002850279102955im

julia> induce_velocity(point, srcs[1])
0.08722212007570912 + 0.14002850279102955im

julia> induce_velocity(z, srcs)
-0.4453372874427177 - 0.10592646656959151im

julia> induce_velocity(point, srcs)
-0.4453372874427177 - 0.10592646656959151im
```
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

# Example

```jldoctest
julia> cluster₁ = Vortex.Point.(rand(Complex128, 5), rand(5));

julia> cluster₂ = Vortex.Point.(rand(Complex128, 5), rand(5));

julia> targets = (cluster₁, cluster₂);

julia> sources = Vortex.Blob.(rand(Complex128), rand(10), 0.1);

julia> ẋs = allocate_velocity(targets);

julia> induce_velocity!(ẋs, targets, sources);

julia> ẋs
(Complex{Float64}[-1.28772-1.82158im, 1.9386-1.64147im, -1.56438+1.57158im, -0.626254+0.375842im, -0.806568-0.213201im], Complex{Float64}[-0.583672-2.26031im, -0.329778-1.43388im, 0.426927+1.55352im, -0.93755+0.241361im, -1.08949-0.35598im])
```
"""
function induce_velocity!(ws::AbstractArray, targets::AbstractArray, source)
    for i in 1:length(targets)
        ws[i] += induce_velocity(targets[i], source)
    end
    ws
end

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

# Example

```jldoctest
julia> points = Vortex.Point.([-1, 1], 1.0)
2-element Array{VortexModel.Vortex.Points.Point,1}:
 Point Vortex: z = -1.0 + 0.0im, Γ = 1.0
 Point Vortex: z = 1.0 + 0.0im, Γ = 1.0

julia> vels = allocate_velocity(points)
2-element Array{Complex{Float64},1}:
 0.0+0.0im
 0.0+0.0im

julia> self_induce_velocity!(vels, points)

julia> vels # should be ±0.25im/π
2-element Array{Complex{Float64},1}:
 0.0-0.0795775im
 0.0+0.0795775im
```
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
    advect(src::PointSource, velocity::Complex128, Δt)

Return a new vortex element that represents `src` advected by `velocity` over `Δt`
If this method is implemented by any type `T <: PointSource`,
then an array of type `AbstractArray{T}` can be passed in the first two arguments of [`advect!`](@ref).

# Example

```jldoctest
julia> point = Vortex.Point(1.0 + 0.0, 1.0);

julia> advect(point, 1.0im, 1e-2)
Point Vortex: z = 1.0 + 0.01im, Γ = 1.0
```
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

# Example

```jldoctest
julia> points₋ = [Vortex.Point(x + 0im, 1.0) for x in 1:5];

julia> points₊ = Vector{Vortex.Point}(5);

julia> vels = [ y*im for y in 1.0:5 ];

julia> advect!(points₊, points₋, vels, 1e-2)

julia> points₊
5-element Array{VortexModel.Vortex.Points.Point,1}:
 Point Vortex: z = 1.0 + 0.01im, Γ = 1.0
 Point Vortex: z = 2.0 + 0.02im, Γ = 1.0
 Point Vortex: z = 3.0 + 0.03im, Γ = 1.0
 Point Vortex: z = 4.0 + 0.04im, Γ = 1.0
 Point Vortex: z = 5.0 + 0.05im, Γ = 1.0
```
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
