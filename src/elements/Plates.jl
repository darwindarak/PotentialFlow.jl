module Plates
using DocStringExtensions

export Plate, bound_circulation, bound_circulation!,
    enforce_no_flow_through!, vorticity_flux, suction_parameters, unit_impulse

import ..Vortex
import ..Vortex:@get, MappedVector
import Base: length, deserialize, AbstractSerializer

include("plates/chebyshev.jl")

"""
    Vortex.Plate <: VortexCompositeSource

An infinitely thin, flat plate, represented as a bound vortex sheet

# Fields
$(FIELDS)

# Constructors
- `Plate(N, L, c, α)`
"""
mutable struct Plate <: Vortex.Element
    "chord length"
    L::Float64
    "centroid"
    c::Complex128
    "centroid velocity"
    α::Float64

    "translational motion parameters"
     ċparams::Vortex.MotionParams

    "rotational motion parameters"
    α̇params::Vortex.MotionParams

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

    "Preplanned discrete Chebyshev transform"
    dchebt!::Chebyshev.Transform{Complex128, true}
end

function deserialize(s::AbstractSerializer, t::Type{Plate})
    fields = Tuple(deserialize(s) for i in 1:11)
    deserialize(s)
    dchebt! = Chebyshev.plan_transform!(fields[9])
    Plate(fields..., dchebt!)
end

function Plate(N, L, c, α, ċparams, α̇params)
    ss = Chebyshev.nodes(N)
    zs = c + 0.5L*ss*exp(im*α)

    C  = zeros(Complex128, N)
    A = MappedVector(imag, C, 1)

    dchebt! = Chebyshev.plan_transform!(C)

    Plate(L, c, α,  ċparams, α̇params, 0.0, N, ss, zs, A, C, 0.0, 0.0, dchebt!)
end

length(p::Plate) = p.N
Vortex.circulation(p::Plate) = p.Γ

# For now, the velocity of the plate is explicitly set
Vortex.allocate_velocity(::Plate) = PlateMotion(0.0, 0.0)
Vortex.self_induce_velocity!(_, ::Plate) = nothing

# Add some motion types
struct Sinusoidal <: Vortex.MotionParams
    A::Float64
    f::Float64
    phi::Float64
    mean::Float64
end

struct SteadyMotion <: Vortex.MotionParams
    v::Float64
end

function motion_function(m::SteadyMotion,t)
   m.v
end

function motion_function(m::Sinusoidal,t)
    @get m (A, f, phi, mean)
    mean + A*sin(2*π*f*t + phi)
end

function motion_function(p::Plate,t)
    PlateMotion(motion_function(p.ċparams,t),
                motion_function(p.α̇params,t))
end

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

function Vortex.induce_velocity!(ws::Vector, p::Plate, sources::T) where T <: Union{Tuple, AbstractArray}
    for source in sources
        Vortex.induce_velocity!(ws, p, source)
    end
    ws
end

@Vortex.kind Plate Vortex.Singleton

function Vortex.induce_velocity(p::Plate, src)
    out = Vortex.allocate_velocity(p.zs)
    Vortex.induce_velocity!(out, p, src)
end

function Vortex.induce_velocity!(ws::Vector, p::Plate, src)
    _singular_velocity!(ws, p, Vortex.unwrap(src),
                        Vortex.kind(Vortex.unwrap_src(src)))
end

function _singular_velocity!(ws, p, src, ::Type{Vortex.Singleton})
    Vortex.induce_velocity!(ws, p.zs, Vortex.Point(src))
end

function _singular_velocity!(ws, p, src, ::Type{Vortex.Group})
    for i in eachindex(src)
        Vortex.induce_velocity!(ws, p, src[i])
    end
    ws
end

mutable struct PlateMotion
    ċ::Complex128
    α̇::Float64
end

Vortex.induce_velocity!(::PlateMotion, target::Plate, source) = nothing
Vortex.reset_velocity!(::PlateMotion, src) = nothing

function Vortex.advect!(plate₊::Plate, plate₋::Plate, ṗ::PlateMotion, Δt)
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
unit_impulse(src, plate::Plate) = unit_impulse(Vortex.position(src), plate)

include("plates/boundary_conditions.jl")
include("plates/circulation.jl")

doc"""
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

    Γs₊ = Vortex.circulation(ambient_sys) .+ bound_circulation(plate)
    Δp .+= (Γs₊ .- Γs₋)./Δt

    Δp, Γs₊
end

function Base.show(io::IO, p::Plate)
    lesp, tesp = suction_parameters(p)
    println(io, "Plate: N = $(p.N), L = $(p.L), c = $(p.c), α = $(round(rad2deg(p.α),2))ᵒ")
    print(io, "       LESP = $(round(lesp,2)), TESP = $(round(tesp,2))")
end

end
