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

position(s::PointSource) = error("`Vortex.position` must be defined for $(typeof(s))")

circulation(vs::Collection) = mapreduce(circulation, +, 0.0, vs)
impulse(vs::Collection) = mapreduce(impulse, +, 0.0, vs)


function allocate_velocity(array::AbstractArray{T}) where {T <: Union{PointSource, Complex128}}
    zeros(Complex128, length(array))
end

allocate_velocity(group::Tuple) = map(allocate_velocity, group)

reset_velocity!(array::AbstractArray{Complex128}) = fill!(array, zero(Complex128))
reset_velocity!(::Void, src) = nothing

function reset_velocity!(group::Tuple)
    foreach(reset_velocity!, group)
    group
end

function reset_velocity!(array::AbstractArray{Complex128}, v)
    resize!(array, length(v))
    fill!(array, zero(Complex128))
end

function reset_velocity!(group::Tuple, src)
    foreach(reset_velocity!, group, src)
    group
end

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

function induce_velocity!(ws::AbstractArray, targets::AbstractArray, source)
    for i in 1:length(targets)
        ws[i] += induce_velocity(targets[i], source)
    end
    ws
end

induce_velocity!(::Void, targets, source) = nothing

function induce_velocity!(ws::Tuple, targets::Tuple, source)
    for i in 1:length(targets)
        induce_velocity!(ws[i], targets[i], source)
    end
    ws
end

function mutually_induce_velocity!(ws₁, ws₂, v₁, v₂)
    induce_velocity!(ws₁, v₁, v₂)
    induce_velocity!(ws₂, v₂, v₁)
    nothing
end

self_induce_velocity!(ws, v) = error("`Vortex.self_induce_velocity!` not defined for $(typeof(v))")

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

advect(t::PointSource, w, Δt) = error("`advect` not defined for $(typeof(v))")
function advect!{T <: Collection}(vs₊::T, vs₋::T, w, Δt)
    for (i, v) in enumerate(vs₋)
        vs₊[i] = advect(v, w[i], Δt)
    end
    nothing
end

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
