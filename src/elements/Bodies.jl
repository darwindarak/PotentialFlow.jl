module Bodies

using DocStringExtensions

export PowerBody, RigidBodyMotions, streamfunction

using ..Points
using ..Blobs

using ..Elements
import ..Elements: position, impulse, circulation
import ..Motions: induce_velocity, induce_velocity!, mutually_induce_velocity!, self_induce_velocity,
                  self_induce_velocity!, allocate_velocity, advect!, reset_velocity!, streamfunction

import ..Utils:@get, MappedVector



include("bodies/RigidBodyMotions.jl")
using .RigidBodyMotions

#=
export bound_circulation, bound_circulation!,
       enforce_no_flow_through!, vorticity_flux, suction_parameters, unit_impulse, force
=#

"""
    PowerBody <: Elements.Element

A two-dimensional body, described by power series mapping

# Constructors
- `PowerBody(C, c, α)`
"""
mutable struct PowerBody <: Element
    "power series coefficients, C[1] -> c₁, C[2] -> c₀, etc"
    C::Vector{Complex128}
    "centroid"
    c::Complex128
    "angle"
    α::Float64

    "number of plotting control points"
    N::Int

    "control point coordinates in circle space"
    ζs::Vector{Complex128}

    "control point coordinates in body-fixed space"
    z̃s::Vector{Complex128}

    "map Jacobian in body-fixed coordinates"
    dz̃dζs::Vector{Complex128}

    "tangent in body-fixed coordinates"
    τs::Vector{Complex128}

    "control point coordinates in inertial space"
    zs::Vector{Complex128}

    "coefficients of power series of |z̃(ζ)|²"
    D::Vector{Complex128}

    "total circulation"
    Γ::Float64

end
@kind PowerBody Singleton


#function deserialize(s::AbstractSerializer, t::Type{Plate})
#    fields = Tuple(deserialize(s) for i in 1:11)
#    deserialize(s)
#    dchebt! = Chebyshev.plan_transform!(fields[9])
#    Plate(fields..., dchebt!)
#end


function PowerBody(C::Vector{Complex128}, c::Complex128, α::Float64)
    N = 200
    ζs = circle(N)
    z̃s = powermap(ζs,C)
    dz̃dζs = d_powermap(ζs,C)
    zs = c + z̃s*exp(im*α)
    τs = zero(ζs)
    @. τs = im*ζs*dz̃dζs/abs(dz̃dζs)

    nc = length(C)
    D = zeros(Complex128,nc+1)
    Cbig = [C;zeros(Complex128,nc)]
    for l = 0:nc, j = 1:nc
      D[l+1] += Cbig[j+l]*conj(Cbig[j])
    end
    D[1] *= 0.5

    PowerBody(C, c, α, N, ζs, z̃s, dz̃dζs, τs, zs, D, 0.0)
end

circle(N) = [exp(im*2π*(i-1)/N) for i in 1:N]

function powermap(ζ::Complex128,C::Vector{Complex128})
  ζⁿ = ζ
  z = zero(ζ)
  for c in C
    z += c*ζⁿ
    ζⁿ /= ζ
  end
  z
end

powermap(ζs::Vector{Complex128},C::Vector{Complex128}) = [powermap(ζ,C) for ζ in ζs]


function d_powermap(ζ::Complex128,C::Vector{Complex128})
  dzdζ = C[1]
  ζⁿ = 1/ζ^2
  for n in 1:length(C)-2
    dzdζ -= n*C[n+2]*ζⁿ
    ζⁿ /= ζ
  end
  dzdζ
end

d_powermap(ζs::Vector{Complex128},C::Vector{Complex128}) = [d_powermap(ζ,C) for ζ in ζs]

jacobian(b::PowerBody) = b.dz̃ds

centroid(b::PowerBody) = b.c

angle(b::PowerBody) = b.α

Base.length(b::PowerBody) = b.N
circulation(b::PowerBody) = b.Γ

Elements.conftransform(z::Complex128,b::PowerBody) = powermap(z,b.C)

Elements.image(z::Complex128,b::PowerBody) = 1.0/conj(z)


Elements.conftransform(s::T,b::PowerBody) where T <: Union{Blob,Point} = Elements.conftransform(s.z,b)


Elements.image(s::T,b::PowerBody) where T <: Union{Blob,Point} = Elements.image(s.z,b)


function get_image(sources::T, b::PowerBody) where T <: Union{Tuple, AbstractArray}
    targ = Points.Point[]
    for source in sources
        push!(targ,get_image(source,b))
    end
    targ
end

function get_image(src, b::PowerBody)
    get_image(Elements.unwrap_src(src), b, kind(Elements.unwrap_src(src)))
end

get_image(src, b::PowerBody, ::Type{Singleton}) = get_image(src,b)

function get_image(src, b::PowerBody, ::Type{Group})
  targ = Points.Point[]
  for i in eachindex(src)
      push!(targ,get_image(src[i], b))
  end
  targ
end

function get_image(src::Union{Blob{T},Point{T}}, b::PowerBody) where T <: Complex
    Point{T}(Elements.image(src.z,b),src.S)
end

function get_image(src::Union{Blob{T},Point{T}}, b::PowerBody) where T <: Real
    Point{T}(Elements.image(src.z,b),-src.S)
end

function allocate_velocity(::PowerBody)
    warn("Body kinematics should be initialized manually.  This simply returns a stationary motion")
    RigidBodyMotion(0.0, 0.0)
end

function self_induce_velocity!(motion, ::PowerBody, t)
    motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)
    motion
end

function induce_velocity(ζ::Complex128, b::PowerBody, m::RigidBodyMotion, t)
    @get b (C, D, α, c)
    self_induce_velocity!(m,b,t);

    @get m (ċ,α̇)

    dz̃dζ = d_powermap(ζ,C)

    c̃̇ = ċ*exp(-im*α)

    ζ⁻ˡ = 1/ζ^2
    w̃ = c̃̇*conj(C[1])*ζ⁻ˡ + conj(c̃̇)*(dz̃dζ-C[1])
    for l = 2:length(D)
        w̃ += im*(l-1)*α̇*D[l]*ζ⁻ˡ
        ζ⁻ˡ /= ζ
    end

    w̃/dz̃dζ

end

function streamfunction(ζ::Complex128, b::PowerBody, m::RigidBodyMotion, t)
  @get b (C, D, α, c)
  self_induce_velocity!(m,b,t);

  @get m (ċ,α̇)

  z̃ = powermap(ζ,C)

  c̃̇ = ċ*exp(-im*α)

  ζ⁻ˡ = 1/ζ
  F = -c̃̇*conj(C[1])*ζ⁻ˡ + conj(c̃̇)*(z̃-C[1]*ζ-C[2]) - im*α̇*D[1]
  for l = 2:length(D)
      F -= im*α̇*D[l]*ζ⁻ˡ
      ζ⁻ˡ /= ζ
  end

  imag(F)

end

function streamfunction(ζ::Complex128, b::PowerBody, src, t)

  srcimg = get_image(src,b)
  sys = (src,srcimg)

  ψ = streamfunction(ζ,sys)

  ψ

end

function streamfunction(ζ::Complex128, b::PowerBody, Winf::Complex128, t)
  @get b (C, D, α, c)

  W̃inf = Winf*exp(im*α)

  F = W̃inf*C[1]*ζ + conj(W̃inf*C[1])/ζ

  imag(F)

end



#=

function impulse(p::Plate)
    @get p (c, B₀, α, Γ, L, A)
    -im*c*Γ - exp(im*α)*π*(0.5L)^2*im*(A[0] - 0.5A[2] - B₀)
end




normal(z, α) = imag(exp(-im*α)*z)
tangent(z, α) = real(exp(-im*α)*z)



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
function reset_velocity!(m::RigidBodyMotion, src)
    m.ċ = m.c̈ = zero(Complex128)
    m.α̇ = zero(Complex128)
    m
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
`src` can be either a `Complex128` or a subtype of `Vortex.PointSource`.
"""
function unit_impulse(z::Complex128, plate::Plate)
    z̃ = 2(z - plate.c)*exp(-im*plate.α)/plate.L
    unit_impulse(z̃)
end
unit_impulse(z̃) = -im*(z̃ + real(√(z̃ - 1)*√(z̃ + 1) - z̃))
unit_impulse(src, plate::Plate) = unit_impulse(Elements.position(src), plate)

include("plates/boundary_conditions.jl")
include("plates/circulation.jl")
include("plates/force.jl")

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
=#

function Base.show(io::IO, b::PowerBody)
    println(io, "Power series body: C = $(b.C), c = $(b.c), α = $(round(rad2deg(b.α),2))ᵒ")
end


end
