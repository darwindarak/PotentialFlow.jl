module Sheets

export Sheet
import ..Vortex
import Base: length

"""
    Vortex.Sheet <: Vortex.CompositeSource

A vortex sheet represented by vortex blob control points

## Fields
- `blobs`: the underlying array of vortex blobs
- `Γs`: the cumulated sum of circulation starting from the first control point
- `δ`: the blob radius of all the vortex blobs

## Constructors:

- `Sheet(zs, Γs, δ)` where `zs` is an array of positions for the control points
"""
mutable struct Sheet <: Vortex.CompositeSource
    blobs::Vector{Vortex.Blob}
    Γs::Vector{Float64}
    δ::Float64
end

length(s::Sheet) = length(s.blobs)
Vortex.circulation(s::Sheet) = s.Γs[end] - s.Γs[1]
Vortex.impulse(s::Sheet) = Vortex.impulse(s.blobs)

Vortex.allocate_velocity(s::Sheet) = zeros(Complex128, length(s.blobs))

function Sheet(zs::AbstractArray{T}, Γs::AbstractArray{Float64}, δ::Float64) where {T <: Number}
    dΓs = compute_trapezoidal_weights(Γs)
    blobs = Vortex.Blob.(zs, dΓs, δ)

    return Sheet(blobs, Γs, δ)
end

for T in Vortex.TargetTypes
    @eval Vortex.induce_velocity(t::$T, s::Sheet) = Vortex.induce_velocity(t, s.blobs)
end
Vortex.induce_velocity(t::Sheet, source) = Vortex.induce_velocity(t.blobs, source)
Vortex.induce_velocity(t::Sheet, s::Sheet) = Vortex.induce_velocity(t.blobs, s.blobs)

Vortex.induce_velocity!(ws::Vector, s::Sheet, source) = Vortex.induce_velocity!(ws, s.blobs, source)

function Vortex.mutually_induce_velocity!(ws₁, ws₂, sheet₁::Sheet, sheet₂::Sheet)
    Vortex.mutually_induce_velocity!(ws₁, ws₂, sheet₁.blobs, sheet₂.blobs)
    nothing
end

function Vortex.self_induce_velocity!(ws, sheet::Sheet)
    Vortex.self_induce_velocity!(ws, sheet.blobs)
end

function compute_trapezoidal_weights(Γs)
    N = length(Γs)

    dΓs = similar(Γs)
    dΓs[1] = 0.5*(Γs[2] - Γs[1])
    for i in 2:N-1
        dΓs[i] = 0.5*(Γs[i+1] - Γs[i-1])
    end
    dΓs[N] = 0.5*(Γs[N] - Γs[N-1])

    return dΓs
end

function Vortex.advect!(sheet₊::Sheet, sheet₋::Sheet, ws, Δt)
    if sheet₊ != sheet₋
        N₊ = length(sheet₊)
        N₋ = length(sheet₋)
        if N₊ != N₋
            resize!(sheet₊.blobs, N₋)
            resize!(sheet₊.Γs, N₋)
        end
        copy!(sheet₊.Γs, sheet₋.Γs)
        sheet₊.δ = sheet₋.δ
    end

    Vortex.advect!(sheet₊.blobs, sheet₋.blobs, ws, Δt)
    nothing
end

function Base.show(io::IO, s::Sheet)
    L = sum(abs, diff(getfield.(s.blobs, :z)))
    print(io, "Vortex Sheet: L ≈ $(round(L, 3)), Γ = $(round(s.Γs[end] - s.Γs[1], 3)), δ = $(round(s.δ, 3))")
end

include("sheets/surgery.jl")

end
