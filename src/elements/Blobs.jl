module Blobs

export Blob

import ..Vortex

"""
An immutable structure representing a vortex blob

## Fields
- `z`: position
- `Γ`: circulation
- `δ`: blob radius
"""
struct Blob <: Vortex.PointSource
    z::Complex128
    Γ::Float64
    δ::Float64
end

Vortex.position(b::Blob) = b.z
Vortex.circulation(b::Blob) = b.Γ
Vortex.impulse(b::Blob) = -im*b.z*b.Γ

blob_kernel(z, δ) = 0.5im*z/(π*(abs2(z) + δ^2))

function Vortex.induce_velocity(z::Complex128, b::Blob)
    b.Γ*blob_kernel(z - b.z, b.δ)
end

function Vortex.induce_velocity(target::Blob, source::Blob)
    δ = √(0.5(target.δ^2 + source.δ^2))
    b.Γ*blob_kernel(target.z - source.z, δ)
end

function Vortex.mutually_induce_velocity!(ws₁, ws₂,
                                          blobs₁::Vector{Blob},
                                          blobs₂::Vector{Blob})
    for (s, source) in enumerate(blobs₁)
        for (t, target) in enumerate(blobs₂)
            δ = √(0.5(target.δ^2 + source.δ^2))
            K = blob_kernel(target.z - source.z, δ)
            ws₂[t] += source.Γ*K
            ws₁[s] -= target.Γ*K
        end
    end
    nothing
end

function Vortex.self_induce_velocity!(ws, blobs::Vector{Blob})
    N = length(blobs)

    for s in 1:N, t in s+1:N
        δ = √(0.5(blobs[t].δ^2 + blobs[s].δ^2))
        K = blob_kernel(blobs[t].z -blobs[s].z, δ)
        ws[t] += blobs[s].Γ*K
        ws[s] -= blobs[t].Γ*K
    end
    nothing
end

Vortex.advect(b::Blob, w::Complex128, Δt::Float64) = Blob(b.z + w*Δt, b.Γ, b.δ)

function Base.show(io::IO, b::Blob)
    print(io, "Vortex Blob: z = $(round(b.z, 3)), Γ = $(round(b.Γ, 3)), δ = $(round(b.δ, 3))")
end

end
