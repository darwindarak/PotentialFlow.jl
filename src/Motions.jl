module Motions

export induce_velocity, induce_velocity!, allocate_velocity,
       self_induce_velocity, self_induce_velocity!,
       mutually_induce_velocity!, reset_velocity!,
       advect, advect!

using ..Properties
using ..Elements
using ..RigidBodyMotions


"""
    induce_velocity(target, element, time)

Compute the velocity induced by `element` on `target`

`target` can be:

- a `ComplexF64`
- a subtype of `Vortex.PointSource`
- an array or tuple of vortex elements

while the `element` can be:

- any subtype of `Vortex.Element`
- an array or tuple of vortex elements

# Example
```jldoctest
julia> z = rand(ComplexF64)
0.23603334566204692 + 0.34651701419196046im

julia> point = Vortex.Point(z, rand());

julia> srcs = Vortex.Point.(rand(ComplexF64, 10), rand(10));

julia> induce_velocity(z, srcs[1], 0.0)
0.08722212007570912 + 0.14002850279102955im

julia> induce_velocity(point, srcs[1], 0.0)
0.08722212007570912 + 0.14002850279102955im

julia> induce_velocity(z, srcs, 0.0)
-0.4453372874427177 - 0.10592646656959151im

julia> induce_velocity(point, srcs, 0.0)
-0.4453372874427177 - 0.10592646656959151im
```
"""
@property begin
    signature = induce_velocity(targ::Target, src::Source, t)
    preallocator = allocate_velocity
    stype = ComplexF64
end

@doc """
    allocate_velocity(srcs)

Allocate arrays of `ComplexF64` to match the structure of `srcs`

# Example

```jldoctest
julia> points = Vortex.Point.(rand(ComplexF64, 2), rand(2));

julia> blobs  = Vortex.Blob.(rand(ComplexF64, 3), rand(3), rand(3));

julia> allocate_velocity(points)
2-element Array{Complex{Float64},1}:
 0.0 + 0.0im
 0.0 + 0.0im

julia> allocate_velocity((points, blobs))
(Complex{Float64}[0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im])
```
""" allocate_velocity

@doc """
    induce_velocity!(vels, target, element, time)

Compute the velocity induced by `element` on `target` and store the result in `vels`

`vels` should be the output of a call to [`allocate_velocity`](@ref),
`target` can be an array or tuple of vortex elements,
while the `element` can be:
- any subtype of `Vortex.Element`
- an array or tuple of vortex elements

# Example

```jldoctest
julia> cluster₁ = Vortex.Point.(rand(ComplexF64, 5), rand(5));

julia> cluster₂ = Vortex.Point.(rand(ComplexF64, 5), rand(5));

julia> targets = (cluster₁, cluster₂);

julia> sources = Vortex.Blob.(rand(ComplexF64), rand(10), 0.1);

julia> ẋs = allocate_velocity(targets);

julia> induce_velocity!(ẋs, targets, sources, 0.0)
(Complex{Float64}[-1.28772-1.82158im, 1.9386-1.64147im, -1.56438+1.57158im, -0.626254+0.375842im, -0.806568-0.213201im], Complex{Float64}[-0.583672-2.26031im, -0.329778-1.43388im, 0.426927+1.55352im, -0.93755+0.241361im, -1.08949-0.35598im])
```
""" induce_velocity!

"""
    mutually_induce_velocity!(vs₁, vs₂, e₁, e₂, t)

Compute the mutually induced velocities between `e₁` and `e₂` at time
`t` and store the results in `vs₁` and `vs₂`

The default implementation simply calls [`induce_velocity!`](@ref)
twice.  This method is meant to be overwritten to take advantage of
symmetries in certain pairwise vortex interations.  For example, the
velocity kernel for a point vortex is antisymmetric, so in computing
the mutually induced velocities of two arrays of point vortices, we
can half the number of calls to the velocity kernel.
"""
function mutually_induce_velocity!(ws₁, ws₂, v₁, v₂, t)
    induce_velocity!(ws₁, v₁, v₂, t)
    induce_velocity!(ws₂, v₂, v₁, t)
    nothing
end

"""
    self_induce_velocity!(vels, elements, time)

Compute the self induced velocity of one or more vortex elements

This involves a recursive call to `self_induce_velocity!` and pairwise
calls to [`mutually_induce_velocity!`](@ref).

# Example

```jldoctest
julia> points = Vortex.Point.([-1, 1], 1.0)
2-element Array{PotentialFlow.Points.Point{Float64},1}:
 Vortex.Point(-1.0 + 0.0im, 1.0)
 Vortex.Point(1.0 + 0.0im, 1.0)

julia> vels = allocate_velocity(points)
2-element Array{Complex{Float64},1}:
 0.0 + 0.0im
 0.0 + 0.0im

julia> self_induce_velocity!(vels, points, 0.0) # should be ±0.25im/π
2-element Array{Complex{Float64},1}:
 0.0 - 0.07957747154594767im
 0.0 + 0.07957747154594767im
```
"""
function self_induce_velocity!(out, group::Tuple, t)
    N = length(group)
    for s in 1:N
        self_induce_velocity!(out[s], group[s], t)
        for t in s+1:N
            mutually_induce_velocity!(out[s], out[t], group[s], group[t], t)
        end
    end
    out
end

function self_induce_velocity!(out, src, t)
    _self_induce_velocity!(out, Elements.unwrap_src(src), kind(Elements.unwrap_src(src)), t)
end

function _self_induce_velocity!(out, src, ::Type{Group}, t)
    N = length(src)
    for s in 1:N
        out[s] += self_induce_velocity(src[s], t)
        for t in s+1:N
            out[s] += induce_velocity(src[s], src[t], t)
            out[t] += induce_velocity(src[t], src[s], t)
        end
    end
    out
end

self_induce_velocity(src, t) = _self_induce_velocity(src, t, kind(Elements.unwrap_src(src)))
_self_induce_velocity(src, t, ::Type{Singleton}) = zero(ComplexF64)
function _self_induce_velocity(src, t, ::Type{Group})
    out = allocate_velocity(src)
    self_induce_velocity!(out, src, t)
end

"""
    reset_velocity!(vels[, srcs])

Set all velocities in `vels` to zero

If `srcs` is provided, then the arrays in `vels` are resized their source counterpart, if necessary.

# Example

```jldoctest
julia> ẋs = (rand(ComplexF64, 1), rand(ComplexF64, 1))
(Complex{Float64}[0.236033+0.346517im], Complex{Float64}[0.312707+0.00790928im])

julia> points = Vortex.Point.(rand(ComplexF64, 2), rand(2));

julia> blobs  = Vortex.Blob.(rand(ComplexF64, 3), rand(3), rand(3));

julia> reset_velocity!(ẋs, (points, blobs));

julia> ẋs
(Complex{Float64}[0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im])
```
"""
reset_velocity!(array::Vector{ComplexF64}) = fill!(array, zero(ComplexF64))

function reset_velocity!(group::Tuple)
    for ẋ in group
        reset_velocity!(ẋ)
    end
    group
end

function reset_velocity!(array::Vector{ComplexF64}, v)
    if length(array) != length(v)
        resize!(array, length(v))
    end
    fill!(array, zero(ComplexF64))
end

function reset_velocity!(group::T, src) where {T <: Tuple}
    for i in 1:length(group)
        reset_velocity!(group[i], src[i])
    end
    group
end

function reset_velocity!(m::RigidBodyMotion, src)
    m.ċ = m.c̈ = zero(ComplexF64)
    m.α̇ = zero(ComplexF64)
    m
end

function reset_velocity!(ẋ::ComplexF64,src)
    ẋ = zero(ComplexF64)
end

"""
    advect(src::Element, velocity::ComplexF64, Δt)

Return a new element that represents `src` advected by `velocity` over
`Δt`.

 If this method is implemented by any type `T` where `kind(T)` is a
`Singleton`, then an array of type `AbstractArray{T}` can be passed in
the first two arguments of [`advect!`](@ref).
Note that this method is usually only defined for singleton elements

# Example

```jldoctest
julia> point = Vortex.Point(1.0 + 0.0, 1.0);

julia> advect(point, 1.0im, 1e-2)
Vortex.Point(1.0 + 0.01im, 1.0)
```
"""
function advect end

"""
    advect!(srcs₊, srcs₋, vels, Δt)

Moves the elements in `srcs₋` by their corresponding velocity in
`vels` over the interval `Δt` and store the results in `src₊`.

`srcs₋` and `srcs₊` can be either a array of vortex elements or a tuple.

# Example

```jldoctest
julia> points₋ = [Vortex.Point(x + 0im, 1.0) for x in 1:5];

julia> points₊ = Vector{Vortex.Point}(undef, 5);

julia> vels = [ y*im for y in 1.0:5 ];

julia> advect!(points₊, points₋, vels, 1e-2);

julia> points₊
5-element Array{PotentialFlow.Points.Point{Float64},1}:
 Vortex.Point(1.0 + 0.01im, 1.0)
 Vortex.Point(2.0 + 0.02im, 1.0)
 Vortex.Point(3.0 + 0.03im, 1.0)
 Vortex.Point(4.0 + 0.04im, 1.0)
 Vortex.Point(5.0 + 0.05im, 1.0)
```
"""
function advect!(vs₊::T, vs₋::T, w, Δt) where
    {T <: Union{AbstractArray, Tuple}}

    for (i, v) in enumerate(vs₋)
        vs₊[i] = advect(v, w[i], Δt)
    end
    vs₊
end

function advect!(group₊::T, group₋::T, ws, Δt) where {T <: Tuple}
    for i in 1:length(group₋)
        advect!(group₊[i], group₋[i], ws[i], Δt)
    end
    group₊
end

end
