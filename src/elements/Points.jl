module Points

export Point

using ..Elements
import ..Motions: induce_velocity, mutually_induce_velocity!, self_induce_velocity!, advect

#== Type definition ==#

"""
    Point <: Elements.Element

An immutable structure representing a point source/vortex

## Fields

- `z`: position
- `S`: strength/circulation
"""
struct Point{T <: Number} <: Element
    z::Complex128
    S::T
    Point{T}(z, s::Real) where T <: Complex = new(z, im*s)
    Point{T}(z, s)       where T <: Complex = new(z, s)
    Point{T}(z, s::Real) where T <: Real = new(z, s)
end
Elements.kind(::Point) = Singleton
Elements.kind(::Type{Point{T}}) where T = Singleton

(p::Point{T})(; z = p.z, S = p.S) where T = Point{T}(z, S)

#== Methods to be extended ==#

Elements.position(p::Point) = p.z

cauchy_kernel(z) = z != zero(z) ? 0.5im/(π*conj(z)) : zero(z)

function induce_velocity(z::Complex128, p::Point, t)
    p.S'*cauchy_kernel(z - p.z)
end

function mutually_induce_velocity!(ws₁, ws₂,
                                   points₁::Vector{Point},
                                   points₂::Vector{Point}, t)
    for (s, source) in enumerate(points₁)
        for (t, target) in enumerate(points₂)
            K = cauchy_kernel(target.z - source.z)
            ws₂[t] += source.Γ*K
            ws₁[s] -= target.Γ*K
        end
    end
    nothing
end

function self_induce_velocity!(ws, points::Vector{Point}, t)
    N = length(points)

    for s in 1:N, t in s+1:N
        K = cauchy_kernel(points[t].z - points[s].z)
        ws[t] += points[s].Γ*K
        ws[s] -= points[t].Γ*K
    end
    ws
end

function advect(p::Point{T}, w::Complex128, Δt::Float64) where T
    Point{T}(p.z + w*Δt, p.S)
end

#function Base.show(io::IO, p::Point)
#    print(io, "Vortex.Point(z = $(round(p.z, 3)), Γ = $(round(p.Γ, 3)))")
#end

end
