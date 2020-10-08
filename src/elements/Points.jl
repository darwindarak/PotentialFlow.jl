module Points

export Point

using ..Elements
import ..Motions: induce_velocity, mutually_induce_velocity!, self_induce_velocity!, advect

#== Type definition ==#

"""
    Point{T,R} <: Elements.Element

An immutable structure representing a point source/vortex

## Fields

- `z::ComplexF64`: position
- `S::T`: strength/circulation
"""
struct Point{T <: Number, R <: Number} <: Element
    z::Complex{R}
    S::T
    Point{T,R}(z, s::Real) where {T <: Complex,R <: Number} = new(z, im*s)
    Point{T,R}(z, s)       where {T <: Complex,R <: Number} = new(z, s)
    Point{T,R}(z, s::Real) where {T <: Real, R <: Number} = new(z, s)
end
#Point(z::Complex{R},S::T) where {T,R} = Point{T,R}(z,S)

Elements.kind(::Point) = Singleton
Elements.kind(::Type{Point{T,R}}) where {T,R} = Singleton

#== Methods to be extended ==#

Elements.position(p::Point) = p.z
Elements.blobradius(p::Point) = 0.0

Elements.streamfunction(z::ComplexF64, p::Point) = real(-0.5p.S*log(z - p.z)/π)

Elements.complexpotential(z::ComplexF64, p::Point) = -0.5im*p.S*log(z - p.z)/π


cauchy_kernel(z) = z != zero(z) ? 0.5im/(π*conj(z)) : zero(z)

function induce_velocity(z::ComplexF64, p::Point, t)
    p.S'*cauchy_kernel(z - p.z)
end


function mutually_induce_velocity!(ws₁, ws₂,
                                   points₁::Vector{Point{T₁}},
                                   points₂::Vector{Point{T₂}}, t) where {T₁, T₂}
    for (s, source) in enumerate(points₁)
        for (t, target) in enumerate(points₂)
            K = cauchy_kernel(target.z - source.z)
            ws₂[t] += source.S'*K
            ws₁[s] -= target.S'*K
        end
    end
    nothing
end

function self_induce_velocity!(ws, points::Vector{Point{T,R}}, t) where {T,R}
    N = length(points)

    for s in 1:N, t in s+1:N
        K = cauchy_kernel(points[t].z - points[s].z)
        ws[t] += points[s].S'*K
        ws[s] -= points[t].S'*K
    end
    ws
end

function advect(p::Point{T,R}, w::ComplexF64, Δt::Float64) where {T,R}
    Point{T,R}(p.z + w*Δt, p.S)
end

#function Base.show(io::IO, p::Point)
#    print(io, "Vortex.Point(z = $(round(p.z, digits=3)), Γ = $(round(p.Γ, digits=3)))")
#end

end
