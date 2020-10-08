module Blobs

export Blob

using ..Elements
import ..Motions: induce_velocity, mutually_induce_velocity!, self_induce_velocity!, advect

#== Type definition ==#

"""
    Blob <: Elements.Element

An immutable structure representing a regularized point source/vortex

## Fields
- `z`: position
- `S`: strength/circulation
- `δ`: regularization radius

## Constructors

"""
struct Blob{T <: Number, R <: Number} <: Element
    z::Complex{R}
    S::T
    δ::Float64
    Blob{T,R}(z, s::Real, δ) where {T <: Complex, R <: Number} = new(z, im*s, δ)
    Blob{T,R}(z, s, δ)       where {T <: Complex, R <: Number} = new(z, s, δ)
    Blob{T,R}(z, s::Real, δ) where {T <: Real, R <: Number}    = new(z, s, δ)
end


Elements.kind(::Blob) = Singleton
Elements.kind(::Type{Blob{T,R}}) where {T,R} = Singleton

#(b::Blob{T,R})(; z = b.z, S = b.S, δ = b.δ) where {T,R} = Blob{T,R}(z, S, δ)

#== Methods to be extended ==#

Elements.position(b::Blob) = b.z
Elements.blobradius(b::Blob) = b.δ
Elements.streamfunction(z::ComplexF64, b::Blob) = real(-0.5b.S*log(z - b.z)/π)

Elements.complexpotential(z::ComplexF64, b::Blob) = -0.5im*b.S*log(z - b.z)/π


blob_kernel(z, δ) = 0.5im*z/(π*(abs2(z) + δ^2))

function induce_velocity(z::ComplexF64, b::Blob, t)
    b.S*blob_kernel(z - b.z, b.δ)
end


function induce_velocity(target::Blob, source::Blob, t)
    δ = √(0.5(target.δ^2 + source.δ^2))
    source.S'*blob_kernel(target.z - source.z, δ)
end

function mutually_induce_velocity!(ws₁, ws₂,
                                   blobs₁::Vector{Blob{T₁}},
                                   blobs₂::Vector{Blob{T₂}}, t) where {T₁, T₂}
    for (s, source) in enumerate(blobs₁)
        for (t, target) in enumerate(blobs₂)
            δ = √(0.5(target.δ^2 + source.δ^2))
            K = blob_kernel(target.z - source.z, δ)
            ws₂[t] += source.S'*K
            ws₁[s] -= target.S'*K
        end
    end
    nothing
end

function self_induce_velocity!(ws, blobs::Vector{Blob{T,R}}, t) where {T,R}
    N = length(blobs)

    for s in 1:N, t in s+1:N
        δ = √(0.5(blobs[t].δ^2 + blobs[s].δ^2))
        K = blob_kernel(blobs[t].z -blobs[s].z, δ)
        ws[t] += blobs[s].S'*K
        ws[s] -= blobs[t].S'*K
    end
    ws
end

function advect(b::Blob{T,R}, w::ComplexF64, Δt::Float64) where {T,R}
    Blob{T,R}(b.z + w*Δt, b.S, b.δ)
end
#
#function Base.show(io::IO, b::Blob)
#    print(io, "Vortex Blob: z = $(round(b.z, digits=3)), Γ = $(round(b.Γ, digits=3)), δ = $(round(b.δ, digits=3))")
#end

end
