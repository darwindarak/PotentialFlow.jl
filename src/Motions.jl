module Motions

export induce_velocity, induce_velocity!, allocate_velocity,
       self_induce_velocity, self_induce_velocity!,
       mutually_induce_velocity!, reset_velocity!,
       advect, advect!

using ..Properties
using ..Elements

@property begin
    signature = induce_velocity(targ::Target, src::Source, t)
    preallocator = allocate_velocity
    stype = Complex128
end

function mutually_induce_velocity!(ws₁, ws₂, v₁, v₂, t)
    induce_velocity!(ws₁, v₁, v₂, t)
    induce_velocity!(ws₂, v₂, v₁, t)
    nothing
end

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
_self_induce_velocity(src, t, ::Type{Singleton}) = zero(Complex128)
function _self_induce_velocity(src, t, ::Type{Group})
    out = allocate_velocity(src)
    self_induce_velocity!(out, src, t)
end

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

function advect!(vs₊::T, vs₋::T, w, Δt) where
    {T <: Union{AbstractArray, Tuple}}

    for (i, v) in enumerate(vs₋)
        vs₊[i] = advect(v, w[i], Δt)
    end
    vs₊
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
    group₊
end

end
