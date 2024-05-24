module Blobs

export Blob

using ..Elements
import ..Motions: induce_velocity, mutually_induce_velocity!, self_induce_velocity!, advect,
                  allocate_velocity, allocate_dveldz, allocate_dveldzstar, dinduce_velocity_dz, dinduce_velocity_dzstar

#== Type definition ==#

"""
    Blob <: Elements.Element

An immutable structure representing a regularized point source/vortex

## Fields
- `z`: position
- `S`: strength/circulation
- `δ`: regularization radius
- `period` : (optional) complex spatial period

## Constructors

"""
struct Blob{S <: Number,R <: Real,P} <: Element
    z::Complex{R}
    S::Union{R,Complex{R}}
    δ::Float64
    period::Number
    Blob{S}(z::Complex{R}, s::T, δ, period)           where {S <: Complex, T <: Real, R <: Real} = (U = promote_type(R,T); new{Complex{U},U,Val{period}}(z, im*convert(U,s),δ,period))
    Blob{S}(z::Complex{R}, s::Complex{T}, δ, period)  where {S <: Complex, T <: Real, R <: Real} = (U = promote_type(R,T); new{Complex{U},U,Val{period}}(z, convert(Complex{U},s),δ,period))
    Blob{S}(z::Complex{R}, s::T, δ, period)           where {S <: Real, T <: Real, R <: Real}    = (U = promote_type(R,T); new{U,U,Val{period}}(z, convert(U,s),δ,period))
end
# struct Blob{T <: Number, R <: Number} <: Element
#     z::Complex{R}
#     S::T
#     δ::Float64
#     Blob{T,R}(z, s::Real, δ) where {T <: Complex, R <: Number} = new(z, im*s, δ)
#     Blob{T,R}(z, s, δ)       where {T <: Complex, R <: Number} = new(z, s, δ)
#     Blob{T,R}(z, s::Real, δ) where {T <: Real, R <: Number}    = new(z, s, δ)
# end


Elements.kind(::Blob) = Singleton
Elements.kind(::Type{Blob{T,R,P}}) where {T,R,P} = Singleton

Elements.property_type(::Blob{T,R,P}) where {T,R,P} = R
Elements.property_type(::Type{Blob{T,R,P}}) where {T,R,P} = R

#Elements.promote_property_type(::Blob{T,R}) where {T,R} = promote_type(T,R)
#Elements.promote_property_type(::Type{Blob{T,R}}) where {T,R} = promote_type(T,R)


#(b::Blob{T,R})(; z = b.z, S = b.S, δ = b.δ) where {T,R} = Blob{T,R}(z, S, δ)

#== Methods to be extended ==#

Elements.position(b::Blob) = b.z
Elements.blobradius(b::Blob) = b.δ

Elements.streamfunction(z::Complex{T}, b::Blob{S,R,Val{Inf}}) where {T,S,R} = real(-0.5b.S*log(z - b.z)/π)
Elements.streamfunction(z::Complex{T}, b::Blob{S,R,Val{P}}) where {T,S,R,P} = real(-0.5b.S*log(im*sin(π*(z - b.z)/P))/π)

Elements.complexpotential(z::Complex{T}, b::Blob{S,R,Val{Inf}}) where {T,S,R} = -0.5im*b.S*log(z - b.z)/π
Elements.complexpotential(z::Complex{T}, b::Blob{S,R,Val{P}}) where {T,S,R,P} = -0.5im*b.S*log(im*sin(π*(z - b.z)/P))/π



# This is the conjugate (K_delta)* of the regularized Cauchy kernel
blob_kernel(z, δ) = 0.5im*z/(π*(abs2(z) + δ^2))

# This is the conjugate (K_delta)* of the regularized Cauchy kernel for period vortices
periodic_blob_kernel(z, δ, λ) = (zp = zmod(z*π/λ)*λ/π; blob_kernel(zp, δ)*zcotz(π*conj(zp)/λ))

zcotz(z) = z != zero(z) ? z*cot(z) : one(z)
zmod(z) = -0.5π + mod(real(z-0.5π),π) + im*imag(z)


# This is d(K_delta)*/dz* = (dK_delta/dz)*
dblob_kernel_dzstar(z, δ) = -0.5im*z^2/(π*(abs2(z) + δ^2)^2)

# This is d(K_delta)*/dz = (dK_delta/dz*)*
dblob_kernel_dz(z, δ) = z != zero(z) ? 0.5im*δ^2/(π*(abs2(z) + δ^2)^2) : zero(z)


# ensures that the velocity is of same type as position
allocate_velocity(v::Vector{Blob{T,R}}) where {T,R} = zeros(Complex{R},length(v))
allocate_dveldz(v::Vector{Blob{T,R}}) where {T,R} = zeros(Complex{R},length(v))
allocate_dveldzstar(v::Vector{Blob{T,R}}) where {T,R} = zeros(Complex{R},length(v))


function induce_velocity(z::Complex{T}, b::Blob{S,R,Val{Inf}}, t) where {T,S,R}
    b.S'*blob_kernel(z - b.z, b.δ)
end

function induce_velocity(z::Complex{T}, b::Blob{S,R,Val{P}}, t) where {T,S,R,P}
    b.S'*periodic_blob_kernel(z - b.z, b.δ, P)
end

function induce_velocity(target::Blob{S,R,Val{Inf}}, source::Blob{S,R,Val{Inf}}, t) where {S,R}
    δ = √(0.5(target.δ^2 + source.δ^2))
    source.S'*blob_kernel(target.z - source.z, δ)
end

function induce_velocity(target::Blob{S₁,R₁,Val{P}}, source::Blob{S₂,R₂,Val{P}}, t) where {S₁,S₂,R₁,R₂,P}
    δ = √(0.5(target.δ^2 + source.δ^2))
    source.S'*periodic_blob_kernel(target.z - source.z, δ, P)
end

function dinduce_velocity_dz(z::Complex{T}, b::Blob, t) where {T}
    b.S*conj(dblob_kernel_dzstar(z - b.z,b.δ))
end

function dinduce_velocity_dzstar(z::Complex{T}, b::Blob, t) where {T}
    b.S*conj(dblob_kernel_dz(z - b.z,b.δ))
end

function mutually_induce_velocity!(ws₁, ws₂,
                                   blobs₁::Vector{Blob{T₁,R₁,Val{Inf}}},
                                   blobs₂::Vector{Blob{T₂,R₂,Val{Inf}}}, t) where {T₁, T₂, R₁, R₂}
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

function mutually_induce_velocity!(ws₁, ws₂,
                                  blobs₁::Vector{Blob{T₁,R₁,Val{P}}},
                                  blobs₂::Vector{Blob{T₂,R₂,Val{P}}}, t) where {T₁, T₂, R₁, R₂,P}
    for (s, source) in enumerate(blobs₁)
        for (t, target) in enumerate(blobs₂)
            δ = √(0.5(target.δ^2 + source.δ^2))
            K = periodic_blob_kernel(target.z - source.z, δ, P)
            ws₂[t] += source.S'*K
            ws₁[s] -= target.S'*K
        end
    end
    nothing
end

function self_induce_velocity!(ws, blobs::Vector{Blob{T,R,Val{Inf}}}, t) where {T,R}
    N = length(blobs)

    for s in 1:N, t in s+1:N
        δ = √(0.5(blobs[t].δ^2 + blobs[s].δ^2))
        K = blob_kernel(blobs[t].z -blobs[s].z, δ)
        ws[t] += blobs[s].S'*K
        ws[s] -= blobs[t].S'*K
    end
    ws
end

function self_induce_velocity!(ws, blobs::Vector{Blob{T,R,Val{P}}}, t) where {T,R,P}
    N = length(blobs)

    for s in 1:N, t in s+1:N
        δ = √(0.5(blobs[t].δ^2 + blobs[s].δ^2))
        K = periodic_blob_kernel(blobs[t].z -blobs[s].z, δ, P)
        ws[t] += blobs[s].S'*K
        ws[s] -= blobs[t].S'*K
    end

    ws
end

function advect(b::Blob{T}, w::ComplexF64, Δt::Float64) where {T}
    Blob{T}(b.z + w*Δt, b.S, b.δ, b.period)
end
#
#function Base.show(io::IO, b::Blob)
#    print(io, "Vortex Blob: z = $(round(b.z, digits=3)), Γ = $(round(b.Γ, digits=3)), δ = $(round(b.δ, digits=3))")
#end

end
