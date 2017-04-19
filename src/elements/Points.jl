module Points

export Point
import ..Vortex

struct Point <: Vortex.PointSource
    z::Complex128
    Γ::Float64
end

Vortex.position(p::Point) = p.z
Vortex.circulation(p::Point) = p.Γ
Vortex.impulse(p::Point) = -im*p.z*p.Γ

cauchy_kernel(z) = 0.5im/(π*conj(z))

function Vortex.induce_velocity(z::Complex128, p::Point)
    p.Γ*cauchy_kernel(z - p.z)
end

function Vortex.mutually_induce_velocity!(ws₁, ws₂,
                                          points₁::Vector{Point},
                                          points₂::Vector{Point})
    for (s, source) in enumerate(points₁)
        for (t, target) in enumerate(points₂)
            K = cauchy_kernel(target.z - source.z)
            ws₂[t] += source.Γ*K
            ws₁[s] -= target.Γ*K
        end
    end
    nothing
end

function Vortex.self_induce_velocity!(ws, points::Vector{Point})
    N = length(points)

    for s in 1:N, t in s+1:N
        K = cauchy_kernel(points[t].z - points[s].z)
        ws[t] += points[s].Γ*K
        ws[s] -= points[t].Γ*K
    end
    nothing
end

Vortex.advect(p::Point, w::Complex128, Δt::Float64) = Point(p.z + w*Δt, p.Γ)

function Base.show(io::IO, b::Point)
    println(io, "Point Vortex: z = $(round(p.z, 3)), Γ = $(round(p.Γ, 3))")
end

end
