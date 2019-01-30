module Plates

using DocStringExtensions
using LinearAlgebra: rmul!

export Plate, bound_circulation, bound_circulation!,
       enforce_no_flow_through!, vorticity_flux, suction_parameters, unit_impulse, force,
       moment

using ..Points
using ..Blobs

using ..Elements
using ..RigidBodyMotions

import ..Elements: position, impulse, circulation
import ..Motions: induce_velocity, induce_velocity!, mutually_induce_velocity!, self_induce_velocity,
                  self_induce_velocity!, allocate_velocity, advect!

import ..Utils:@get, MappedVector

include("plates/chebyshev.jl")

"""
    Plate <: Elements.Element

An infinitely thin, flat plate, represented as a bound vortex sheet

# Constructors
- `Plate(N, L, c, α)`
"""
mutable struct Plate <: Element
    "chord length"
    L::Float64
    "centroid"
    c::ComplexF64
    "centroid velocity"
    α::Float64

    "total circulation"
    Γ::Float64

    "number of control points"
    N::Int
    "normalized positions (within [-1, 1]) of the control points"
    ss::Vector{Float64}
    "control point coordinates"
    zs::Vector{ComplexF64}
    "Chebyshev coefficients of the normal component of velocity induced along the plate by ambient vorticity"
    A::MappedVector{Float64, Vector{ComplexF64}, typeof(imag)}
    "Chebyshev coefficients of the velocity induced along the plate by ambient vorticity"
    C::Vector{ComplexF64}
    "zeroth Chebyshev coefficient associated with body motion"
    B₀::Float64
    "first Chebyshev coefficient associated with body motion"
    B₁::Float64

    "Preplanned discrete Chebyshev transform"
    dchebt!::Chebyshev.Transform{ComplexF64, true}
end
@kind Plate Singleton

#function deserialize(s::AbstractSerializer, t::Type{Plate})
#    fields = Tuple(deserialize(s) for i in 1:11)
#    deserialize(s)
#    dchebt! = Chebyshev.plan_transform!(fields[9])
#    Plate(fields..., dchebt!)
#end

function Plate(N, L, c, α)
    ss = Chebyshev.nodes(N)
    zs = c .+ 0.5L*ss*exp(im*α)

    C  = zeros(ComplexF64, N)
    A = MappedVector(imag, C, 1)

    dchebt! = Chebyshev.plan_transform!(C)

    Plate(L, c, α, 0.0, N, ss, zs, A, C, 0.0, 0.0, dchebt!)
end

Base.length(p::Plate) = p.N

circulation(p::Plate) = p.Γ
function impulse(p::Plate)
    @get p (c, B₀, α, Γ, L, A)
    -im*c*Γ - exp(im*α)*π*(0.5L)^2*im*(A[0] - 0.5A[2] - B₀)
end
function angularimpulse(p::Plate)
    @get p (c, B₀, α, Γ, L, A)
    -0.5*(c*conj(c)+L^2/8)*Γ - real(exp(im*α)*conj(c)*π*(0.5L)^2*(A[0] - 0.5A[2] - B₀)) -
          π*(0.25L)^3*(A[1]-A[3]-B₁)
end

function allocate_velocity(::Plate)
    @warn("Plate kinematics should be initialized manually.  This simply returns a stationary motion")
    RigidBodyMotion(0.0, 0.0)
end
function self_induce_velocity!(motion, ::Plate, t)
    motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)
    motion
end

normal(z, α) = imag(exp(-im*α)*z)
tangent(z, α) = real(exp(-im*α)*z)

function induce_velocity(z::ComplexF64, p::Plate, t)
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

      # note that this provides u+iv (JDE)
    0.5im*w*exp(im*α)
end

function induce_velocity!(ws::Vector, p::Plate, sources::T, t) where T <: Union{Tuple, AbstractArray}
    for source in sources
        induce_velocity!(ws, p, source, t)
    end
    ws
end
function induce_velocity(p::Plate, src, t)
    out = allocate_velocity(p.zs)
    induce_velocity!(out, p, src, t)
end

function induce_velocity!(ws::Vector, p::Plate, src, t)
    _singular_velocity!(ws, p, Elements.unwrap(src), t,
                        kind(Elements.unwrap_src(src)))
end

function _singular_velocity!(ws, p, src::Blob{T}, t, ::Type{Singleton}) where T
    induce_velocity!(ws, p.zs, Point{T}(src.z, src.S), t)
end

function _singular_velocity!(ws, p, src, t, ::Type{Singleton})
    induce_velocity!(ws, p.zs, src, t)
end

function _singular_velocity!(ws, p, src, t, ::Type{Group})
    for i in eachindex(src)
        induce_velocity!(ws, p, src[i], t)
    end
    ws
end

induce_velocity!(m::RigidBodyMotion, target::Plate, source, t) = m


function Elements.streamfunction(z::ComplexF64, p::Plate)
    @get p (N, L, c, α, Γ, A, B₀, B₁)
    z̃ = 2*(z - c)*exp(-im*α)/L
    J = z̃ - √(z̃ - 1)*√(z̃ + 1)

    ψ  = -real(J)*2(A[0] - B₀)
    ψ -= (log(abs(J)) + log(2))*(A[1] - B₁ + 2Γ/(π*L))
    ψ -= (real(0.25J^2) - 0.5(log(abs(J)) + log(2)))*(A[1] - B₁)

    for n in 2:N-1
        ψ -= A[n]*0.5real(J^(n+1)/(n+1) - J^(n-1)/(n-1))
    end

    return -ψ/4
end

function advect!(plate₊::Plate, plate₋::Plate, ṗ::RigidBodyMotion, Δt)
    if plate₊ != plate₋
        plate₊.L    = plate₋.L
        plate₊.Γ    = plate₋.Γ
        plate₊.B₀   = plate₋.B₀
        plate₊.B₁   = plate₋.B₁
        if plate₊.N != plate₋.N
            plate₊.N    = plate₋.N
            resize!(plate₊.ss, plate₊.N)
            resize!(plate₊.zs, plate₊.N)
            resize!(plate₊.C, plate₊.N)
        end
        plate₊.ss   .= plate₋.ss
        plate₊.zs   .= plate₋.zs
        plate₊.C    .= plate₋.C
        plate₊.dchebt! = plate₋.dchebt!
    end
    plate₊.c = plate₋.c + ṗ.ċ*Δt
    plate₊.α = plate₋.α + ṗ.α̇*Δt

    @get plate₊ (c, L, α)

    @. plate₊.zs = c + 0.5L*exp(im*α)*plate₊.ss

    plate₊
end

"""
    unit_impulse(src, plate::Plate)

Compute the impulse per unit circulation of `src` and its associated bound vortex sheet on `plate` (its image vortex)
`src` can be either a `ComplexF64` or a subtype of `Vortex.PointSource`.
"""
function unit_impulse(z::ComplexF64, plate::Plate)
    z̃ = 2(z - plate.c)*exp(-im*plate.α)/plate.L
    unit_impulse(z̃)
end
unit_impulse(z̃) = -im*(z̃ + real(√(z̃ - 1)*√(z̃ + 1) - z̃))
unit_impulse(src, plate::Plate) = unit_impulse(Elements.position(src), plate)

include("plates/boundary_conditions.jl")
include("plates/circulation.jl")
include("plates/force.jl")

raw"""
    surface_pressure(plate, motion, te_sys, Γs₋, Δt)

Compute the pressure difference across the plate along Chebyshev nodes.

!!! note
    The pressure difference across the bound vortex sheet is given by:
    ```math
        [p]_-^+
      = -\rho \left[ \frac{1}{2}(\boldsymbol{v}^+ + \boldsymbol{v}^-)
                   - \boldsymbol{v}_b
             \right]
             \cdot ( \boldsymbol{\gamma} \times \boldsymbol{\hat{n}})
        +\rho \frac{\mathrm{d}\Gamma}{\mathrm{d}t}
    ```
    where ``\rho`` is the fluid density, ``\boldsymbol{v}^\pm`` is the
    velocity on either side of the plate, ``\boldsymbol{v}_b`` is the local
    velocity of the plate, ``\boldsymbol{\gamma}`` is the bound vortex
    sheet strength, and ``\Gamma`` is the integrated circulation.
    We will compute ``\frac{\mathrm{d}\Gamma}{\mathrm{d}t}`` using finite
    differences.  So we will need the circulation along the plate from a
    previous time-step in order to compute the current pressure
    distribution.  We assume that value of circulation at the trailing
    edge of the plate is equal the the net circulation of all the vorticity
    that has been shed from the trailing edge.

# Arguments

- `plate`: we assume that the `Plate` structure that is passed in
  already enforces the no-flow-through condition
- `motion`: the motion of the plate used to compute ``\boldsymbol{v}_b``
- `te_sys`: the system of vortex elements representing the vorticity
  shed from the trailing edge of the plate
- `Γs₋`: the circulation along the plate's Chebyshev nodes, this
  should be equivalent to calling
  `Vortex.circulation(te_sys) .+ Vortex.bound_circulation(plate)`
  from a previous time-step.
- `Δt`: time-step used to compute ``\frac{\mathrm{d}\Gamma}{\mathrm{d}t}``
  using finite differences

# Returns

- `Δp`: the pressure difference across the plate along Chebyshev nodes
- `Γs₊`: the circulation along the plate at the current time-step
  (this value is used in computing the current `Δp` and can be used as
  the `Γs₋` for computing pressure differences at the **next** time-step)
"""
function surface_pressure(plate, motion, ambient_sys, Γs₋, Δt)
    @get plate (C, ss, α)

    Δp = strength(plate) .* (Chebyshev.firstkind(real.(C), ss) .- tangent(motion.ċ, α))

    Γs₊ = circulation(ambient_sys) .+ bound_circulation(plate)
    Δp .+= (Γs₊ .- Γs₋)./Δt

    Δp, Γs₊
end

"""
    edges(plate)

Return the coordinates of the leading and trailing edges

# Example

```jldoctest
julia> p = Plate(128, 1.0, 0, π/4)
Plate: N = 128, L = 1.0, c = 0.0 + 0.0im, α = 45.0ᵒ
       LESP = 0.0, TESP = 0.0

julia> Plates.edges(p)
(0.3535533905932738 + 0.35355339059327373im, -0.3535533905932738 - 0.35355339059327373im)
```
"""
edges(plate) = plate.zs[end], plate.zs[1]

include("plates/pressure.jl")

function Base.show(io::IO, p::Plate)
    lesp, tesp = suction_parameters(p)
    println(io, "Plate: N = $(p.N), L = $(p.L), c = $(p.c), α = $(round(rad2deg(p.α); digits=2))ᵒ")
    print(io, "       LESP = $(round(lesp; digits=2)), TESP = $(round(tesp; digits=2))")
end

end
