module Plates

export Plate, bound_circulation, bound_circulation!,
    enforce_no_flow_through!, vorticity_flux, suction_parameters

import ..Vortex
import ..Vortex.@get
import Base: length

mutable struct Plate <: Vortex.CompositeSource
    L::Float64
    c::Complex128
    ċ::Complex128
    α::Float64
    α̇::Float64

    Γ::Float64

    N::Int
    ss::Vector{Float64}
    zs::Vector{Complex128}
    A::Vector{Float64}
    C::Vector{Complex128}
end

function Plate(N, L, c, α, ċ = zero(Complex128), α̇ = 0.0, c̈ = zero(Complex128), α̈ = 0.0)
    ss = Float64[cos(θ) for θ in linspace(π, 0, N)]
    zs = c + 0.5L*ss*exp(im*α)

    A  = zeros(Float64, N)
    C  = zeros(Complex128, N)

    Plate(L, c, ċ, α, α̇, 0.0, N, ss, zs, A, C)
end


length(p::Plate) = p.N
Vortex.circulation(p::Plate) = p.Γ

# For now, the velocity of the plate is explicitly set
Vortex.allocate_velocity(::Plate) = nothing
Vortex.self_induce_velocity!(::Void, ::Plate) = nothing

function Vortex.impulse(p::Plate)
    @get p (c, ċ, α, Γ, L, A)
    B₀ = normal(ċ, α)
    -im*c*Γ - exp(im*α)*π*0.5L*im*(A[1] - 0.5A[2] - B₀)
end

normal(z, α) = imag(exp(-im*α)*z)
tangent(z, α) = real(exp(-im*α)*z)

function Vortex.induce_velocity(z::Complex128, p::Plate)
    @get p (α, L, c, ċ, α̇, Γ, A)

    z̃ = conj(2*(z - c)*exp(-im*α)/L)

    ρ = √(z̃ - 1)*√(z̃ + 1)
    J = z̃ - ρ

    B₀ = normal(ċ, α)
    B₁ = 0.5α̇*L

    w = (A[2] + 2Γ/(π*L))
    w += 2(A[1] - B₀)*J
    w -= B₁*J^2

    w /= ρ

    Jⁿ = J
    for n in 2:length(A)
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

function Vortex.advect!(plate₊::Plate, plate₋::Plate, ws, Δt)

    if plate₊ != plate₋
        plate₊.L  = plate₋.L 
        plate₊.ċ  = plate₋.ċ 
        plate₊.α̇  = plate₋.α̇ 
        plate₊.Γ  = plate₋.Γ 
        plate₊.N  = plate₋.N 
        plate₊.ss = plate₋.ss 
        plate₊.zs = plate₋.zs
        plate₊.A  = plate₋.A 
        plate₊.C  = plate₋.C 
    end
    plate₊.c = plate₋.c + plate₋.ċ*Δt
    plate₊.α = plate₋.α + plate₋.α̇*Δt
    plate₊.zs .+= plate₊.ċ*Δt
    nothing
end

include("plates/chebyshev.jl")
include("plates/boundary_conditions.jl")
include("plates/circulation.jl")

function Base.show(io::IO, p::Plate)
    println(io, "Plate: N = $(p.N), L = $(p.L)")
    println(io, "       c = $(p.c), α = $(round(rad2deg(p.α),2))ᵒ")
    print(io, "       ċ = $(p.ċ), α̇ = $(round(rad2deg(p.α̇),2))ᵒ/T")
end

end
