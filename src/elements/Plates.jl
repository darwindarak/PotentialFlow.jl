module Plates
using DocStringExtensions

export Plate, bound_circulation, bound_circulation!,
    enforce_no_flow_through!, vorticity_flux, suction_parameters, unit_impulse

import ..Vortex
import ..Vortex:@get, MappedVector
import Base: length

"""
An infinitely thin, flat plate, represented as a bound vortex sheet

# Fields
$(FIELDS)

# Constructors
- `Plate(N, L, c, α)
"""
mutable struct Plate <: Vortex.CompositeSource
    "chord length"
    L::Float64
    "centroid"
    c::Complex128
    "centroid velocity"
    α::Float64

    "total circulation"
    Γ::Float64

    "number of control points"
    N::Int
    "normalized positions (within [-1, 1]) of the control points"
    ss::Vector{Float64}
    "control point coordinates"
    zs::Vector{Complex128}
    "Chebyshev coefficients of the normal component of velocity induced along the plate by ambient vorticity"
    A::MappedVector{Float64, Vector{Complex128}, typeof(imag)}
    "Chebyshev coefficients of the velocity induced along the plate by ambient vorticity"
    C::Vector{Complex128}
    "zeroth Chebyshev coefficient associated with body motion"
    B₀::Float64
    "first Chebyshev coefficient associated with body motion"
    B₁::Float64
end

function Plate(N, L, c, α)
    ss = Float64[cos(θ) for θ in linspace(π, 0, N)]
    zs = c + 0.5L*ss*exp(im*α)

    C  = zeros(Complex128, N)
    A = MappedVector(imag, C, Float64, 1)

    Plate(L, c, α, 0.0, N, ss, zs, A, C, 0.0, 0.0)
end

length(p::Plate) = p.N
Vortex.circulation(p::Plate) = p.Γ

# For now, the velocity of the plate is explicitly set
Vortex.allocate_velocity(::Plate) = PlateMotion(0.0, 0.0)
Vortex.self_induce_velocity!(_, ::Plate) = nothing

function Vortex.impulse(p::Plate)
    @get p (c, B₀, α, Γ, L, A)
    -im*c*Γ - exp(im*α)*π*(0.5L)^2*im*(A[0] - 0.5A[2] - B₀)
end

normal(z, α) = imag(exp(-im*α)*z)
tangent(z, α) = real(exp(-im*α)*z)

function Vortex.induce_velocity(z::Complex128, p::Plate)
    @get p (α, L, c, B₀, B₁, Γ, A)

    z̃ = conj(2*(z - c)*exp(-im*α)/L)

    ρ = √(z̃ - 1)*√(z̃ + 1)
    J = z̃ - ρ

    w = (A[1] + 2Γ/(π*L))
    w += 2(A[0] - B₀)*J
    w -= B₁*J^2

    w /= ρ

    Jⁿ = J
    for n in 1:length(A)-1
        w -= 2A[n]*Jⁿ
        Jⁿ *= J
    end

    0.5im*w*exp(im*α)
end

function Vortex.induce_velocity!(ws::Vector, p::Plate, sources::Vortex.Collection)
    for source in sources
        Vortex.induce_velocity!(ws, p, source)
    end
    ws
end

Vortex.induce_velocity!(ws::Vector, p::Plate, ps::Vortex.PointSource) = Vortex.induce_velocity!(ws, p.zs, ps)

function Vortex.induce_velocity!(ws::Vector, p::Plate, b::Vortex.Blob)
    for (i, z) in enumerate(p.zs)
        ws[i] += b.Γ*Vortex.Points.cauchy_kernel(z - b.z)
    end
    ws
end

function Vortex.induce_velocity!(ws::Vector, p::Plate, blobs::Vector{Vortex.Blob})
    for (i, z) in enumerate(p.zs), b in blobs
        ws[i] += b.Γ*Vortex.Points.cauchy_kernel(z - b.z)
    end
    ws
end

function Vortex.induce_velocity!(ws::Vector, p::Plate, s::Vortex.Sheet)
    Vortex.induce_velocity!(ws, p, s.blobs)
end

mutable struct PlateMotion
    ċ::Complex128
    α̇::Float64
end

Vortex.induce_velocity!(::PlateMotion, target::Plate, source) = nothing
Vortex.reset_velocity!(::PlateMotion, src) = nothing

function Vortex.advect!(plate₊::Plate, plate₋::Plate, ṗ::PlateMotion, Δt)

    if plate₊ != plate₋
        plate₊.L  = plate₋.L 
        plate₊.Γ  = plate₋.Γ 
        plate₊.N  = plate₋.N 
        plate₊.ss = plate₋.ss 
        plate₊.zs = plate₋.zs
        plate₊.A  = plate₋.A 
        plate₊.C  = plate₋.C 
        plate₊.B₀ = plate₋.B₀
        plate₊.B₁ = plate₋.B₁ 
    end
    plate₊.c = plate₋.c + ṗ.ċ*Δt
    plate₊.α = plate₋.α + ṗ.α̇*Δt
    plate₊.zs .+= ṗ.ċ*Δt
    nothing
end

"""
    unit_impulse(src, plate::Plate)

Compute the impulse per unit circulation of `src` and its associated bound vortex sheet on `plate` (its image vortex)
`src` can be either a `Complex128` or a subtype of `Vortex.PointSource`.
"""
function unit_impulse(z::Complex128, plate::Plate)
    z̃ = 2(z - plate.c)*exp(-im*plate.α)/plate.L
    unit_impulse(z̃)
end
unit_impulse(z̃) = -im*(z̃ + real(√(z̃ - 1)*√(z̃ + 1) - z̃))
unit_impulse(src::Vortex.PointSource, plate::Plate) = unit_impulse(Vortex.position(src), plate)

include("plates/chebyshev.jl")
include("plates/boundary_conditions.jl")
include("plates/circulation.jl")

function Base.show(io::IO, p::Plate)
    print(io, "Plate: N = $(p.N), L = $(p.L), c = $(p.c), α = $(round(rad2deg(p.α),2))ᵒ")
end

end
