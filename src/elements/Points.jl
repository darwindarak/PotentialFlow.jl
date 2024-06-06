module Points

export Point

using ..Elements
import ..Motions: induce_velocity, mutually_induce_velocity!, self_induce_velocity!, advect,
                  allocate_velocity, allocate_dveldz, allocate_dveldzstar, dinduce_velocity_dz, dinduce_velocity_dzstar

#== Type definition ==#

"""
    Point{T,R} <: Elements.Element

An immutable structure representing a point source/vortex

## Fields

- `z::ComplexF64`: position
- `S::T`: strength/circulation
- `period` : (optional) complex spatial period

"""
struct Point{S <: Number,R <: Real, P} <: Element
    z::Complex{R}
    S::Union{R,Complex{R}}
    period::Number
    Point{S}(z::Complex{R}, s::T, period)           where {S <: Complex, T <: Real, R <: Real} = (U = promote_type(R,T); new{Complex{U},U,Val{period}}(z, im*convert(U,s),period))
    Point{S}(z::Complex{R}, s::Complex{T}, period)  where {S <: Complex, T <: Real, R <: Real} = (U = promote_type(R,T); new{Complex{U},U,Val{period}}(z, convert(Complex{U},s),period))
    Point{S}(z::Complex{R}, s::T, period)           where {S <: Real, T <: Real, R <: Real}    = (U = promote_type(R,T); new{U,U,Val{period}}(z, convert(U,s),period))
end
# struct Point{T <: Number, R <: Number} <: Element
#     z::Complex{R}
#     S::T
#     Point{T,R}(z, s::Real) where {T <: Complex,R <: Number} = new(z, im*s)
#     Point{T,R}(z, s)       where {T <: Complex,R <: Number} = new(z, s)
#     Point{T,R}(z, s::Real) where {T <: Real, R <: Number} = new(z, s)
# end

#Point(z::Complex{R},S::T) where {T,R} = Point{T,R}(z,S)

Elements.kind(::Point) = Singleton
Elements.kind(::Type{Point{T,R,P}}) where {T,R,P} = Singleton

Elements.property_type(::Point{T,R,P}) where {T,R,P} = R
Elements.property_type(::Type{Point{T,R,P}}) where {T,R,P} = R

#Elements.promote_property_type(::Point{T,R}) where {T,R} = promote_type(T,R)
#Elements.promote_property_type(::Type{Point{T,R}}) where {T,R} = promote_type(T,R)

#== Methods to be extended ==#

Elements.position(p::Point) = p.z
Elements.blobradius(p::Point) = 0.0

Elements.streamfunction(z::Complex{T}, p::Point{S,R,Val{Inf}}) where {T,S,R} = real(-0.5p.S*log(z - p.z)/π)
Elements.streamfunction(z::Complex{T}, b::Point{S,R,Val{P}}) where {T,S,R,P} = real(-0.5b.S*log(im*sin(π*(z - b.z)/P))/π)

Elements.complexpotential(z::Complex{T}, p::Point{S,R,Val{Inf}}) where {T,S,R} = -0.5im*p.S*log(z - p.z)/π
Elements.complexpotential(z::Complex{T}, b::Point{S,R,Val{P}}) where {T,S,R,P} = -0.5im*b.S*log(im*sin(π*(z - b.z)/P))/π


# This is actually the conjugate K* of the Cauchy kernel
cauchy_kernel(z) = z != zero(z) ? 0.5im/(π*conj(z)) : zero(z)

# This is the conjugate K* of the Cauchy kernel for period vortices
periodic_cauchy_kernel(z, λ) = mod(abs(z)/λ,1) != zero(z) ? 0.5*im/λ*cot(π*conj(z)/λ) : zero(z)


# This is dK*/dz* = (dK/dz)*
dcauchy_kernel_dzstar(z) = z != zero(z) ? -0.5im/(π*conj(z)^2) : zero(z)

# ensures that the velocity is of same type as position
allocate_velocity(v::Vector{Point{T,R,P}}) where {T,R,P} = zeros(Complex{R},length(v))
allocate_dveldz(v::Vector{Point{T,R,P}}) where {T,R,P} = zeros(Complex{R},length(v))
allocate_dveldzstar(v::Vector{Point{T,R,P}}) where {T,R,P} = zeros(Complex{R},length(v))

function induce_velocity(z::Complex{T}, p::Point{S,R,Val{Inf}}, t) where {T,S,R}
    p.S'*cauchy_kernel(z - p.z)
end

function induce_velocity(z::Complex{T}, p::Point{S,R,Val{P}}, t) where {T,S,R,P}
    p.S'*periodic_cauchy_kernel(z - p.z,P)
end

function dinduce_velocity_dz(z::Complex{T}, p::Point, t) where {T}
    p.S*conj(dcauchy_kernel_dzstar(z - p.z))
end

dinduce_velocity_dzstar(z::Complex{T}, p::Point, t) where {T} = zero(z)

function mutually_induce_velocity!(ws₁, ws₂,
                                   points₁::Vector{Point{T₁,R₁,Val{Inf}}},
                                   points₂::Vector{Point{T₂,R₂,Val{Inf}}}, t) where {T₁, T₂, R₁, R₂}
    for (s, source) in enumerate(points₁)
        for (t, target) in enumerate(points₂)
            K = cauchy_kernel(target.z - source.z)
            ws₂[t] += source.S'*K
            ws₁[s] -= target.S'*K
        end
    end
    nothing
end

function mutually_induce_velocity!(ws₁, ws₂,
                                    points₁::Vector{Point{T₁,R₁,Val{P}}},
                                    points₂::Vector{Point{T₂,R₂,Val{P}}}, t) where {T₁, T₂, R₁, R₂, P}
    for (s, source) in enumerate(points₁)
        for (t, target) in enumerate(points₂)
            K = periodic_cauchy_kernel(target.z - source.z,P)
            ws₂[t] += source.S'*K
            ws₁[s] -= target.S'*K
        end
        end
    nothing
end

function self_induce_velocity!(ws, points::Vector{Point{T,R,Val{Inf}}}, t) where {T,R}
    N = length(points)

    for s in 1:N, t in s+1:N
        K = cauchy_kernel(points[t].z - points[s].z)
        ws[t] += points[s].S'*K
        ws[s] -= points[t].S'*K
    end
    ws
end

function self_induce_velocity!(ws, points::Vector{Point{T,R,Val{P}}}, t) where {T,R,P}
    N = length(points)

    for s in 1:N, t in s+1:N
        K = periodic_cauchy_kernel(points[t].z - points[s].z,P)
        ws[t] += points[s].S'*K
        ws[s] -= points[t].S'*K
    end
    ws
end

function advect(p::Point{T}, w::ComplexF64, Δt::Float64) where {T}
    Point{T}(p.z + w*Δt, p.S, p.period)
end

#function Base.show(io::IO, p::Point)
#    print(io, "Vortex.Point(z = $(round(p.z, digits=3)), Γ = $(round(p.Γ, digits=3)))")
#end

end
