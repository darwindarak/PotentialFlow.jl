module Motions

export Motion

import ForwardDiff
import Base: +, *, -, >>, <<

abstract type Kinematics end
abstract type Profile end

mutable struct Motion
    ċ::Complex128
    c̈::Complex128
    α̇::Float64

    kin::Kinematics
end

Motion(ċ, α̇) = Motion(complex(ċ), 0.0im, float(α̇), Constant(ċ, α̇))
Motion(kin::Kinematics) = Motion(kin(0)..., kin)

function Base.show(io::IO, m::Motion)
    println(io, "Plate Motion:")
    println(io, "  ċ = $(round(m.ċ, 2))")
    println(io, "  c̈ = $(round(m.c̈, 2))")
    println(io, "  α̇ = $(round(m.α̇, 2))")
    print(io, "  $(m.kin)")
end

struct Constant{C <: Complex, A <: Real} <: Kinematics
    ċ::C
    α̇::A
end
Constant(ċ, α̇) = Constant(complex(ċ), α̇)
(c::Constant{C})(t) where C = c.ċ, zero(C), c.α̇
Base.show(io::IO, c::Constant) = print(io, "Constant (ċ = $(c.ċ), α̇ = $(c.α̇))")

struct Pitchup <: Kinematics
    "Freestream velocity"
    U₀::Float64
    "Axis of rotation, relative to the plate centroid"
    a::Float64

    "Non-dimensional pitch rate ``K = \dot{\alpha}_0\frac{c}{2U_0}``"
    K::Float64

    "Initial angle of attack"
    α₀::Float64
    "Nominal start of pitch up"
    t₀::Float64

    "Total pitching angle"
    Δα::Float64

    α::Profile
    α̇::Profile
    α̈::Profile
end

function Pitchup(U₀, a, K, α₀, t₀, Δα, ramp)
    Δt = 0.5Δα/K
    p = 2K*((ramp >> t₀) - (ramp >> (t₀ + Δt)))
    ṗ = d_dt(p)
    p̈ = d_dt(ṗ)

    Pitchup(U₀, a, K, α₀, t₀, Δα, p, ṗ, p̈)
end

function (p::Pitchup)(t)
    α = p.α₀ + p.α(t)
    α̇ = p.α̇(t)
    α̈ = p.α̈(t)

    ċ = p.U₀ - p.a*im*α̇*exp(im*α)
    c̈ = p.a*exp(im*α)*(α̇^2 - im*α̈)

    return ċ, c̈, α̇
end


struct DerivativeProfile{P} <: Profile
    p::P
end

function Base.show(io::IO, ṗ::DerivativeProfile)
    print(io, "d/dt ($(ṗ.p))")
end

(ṗ::DerivativeProfile)(t) = ForwardDiff.derivative(ṗ.p, t)

d_dt(p::Profile) = DerivativeProfile(p)

struct ScaledProfile{N <: Real, P <: Profile} <: Profile
    s::N
    p::P
end
function Base.show(io::IO, p::ScaledProfile)
    print(io, "$(p.s) × $(p.p)")
end
s::Number * p::Profile = ScaledProfile(s, p)
-(p::Profile) = ScaledProfile(-1, p)
(p::ScaledProfile)(t) = p.s*p.p(t)

struct ShiftedProfile{N <: Real, P <: Profile} <: Profile
    Δt::N
    p::P
end
function Base.show(io::IO, p::ShiftedProfile)
    print(io, "$(p.p) >> $(p.Δt)")
end

(p::ShiftedProfile)(t) = p.p(t - p.Δt)
p::Profile >> Δt::Number = ShiftedProfile(Δt, p)
p::Profile << Δt::Number = ShiftedProfile(-Δt, p)

struct AddedProfiles{T <: Tuple} <: Profile
    ps::T
end
function Base.show(io::IO, Σp::AddedProfiles)
    println(io, "AddedProfiles:")
    for p in Σp.ps
        println(io, "  $p")
    end
end

+(p::Profile, Σp::AddedProfiles) = AddedProfiles((Σp..., p))
+(Σp::AddedProfiles, p::Profile) = AddedProfiles((p, Σp...))
function +(Σp₁::AddedProfiles, Σp₂::AddedProfiles)
    AddedProfiles((Σp₁..., Σp₂...))
end

-(p₁::Profile, p₂::Profile) = p₁ + (-p₂)

+(p::Profile...) = AddedProfiles(p)

function (Σp::AddedProfiles)(t)
    f = 0.0
    for p in Σp.ps
        f += p(t)
    end
    f
end

struct EldredgeRamp <: Profile
    aₛ::Float64
end

function (r::EldredgeRamp)(t)
    0.5(log(2cosh(r.aₛ*t)) + r.aₛ*t)/r.aₛ
end

struct ColoniusRamp <: Profile
    n::Int
end

function (r::ColoniusRamp)(t)
    Δt = t + 0.5
    if Δt ≤ 0
        0.0
    elseif Δt ≥ 1
        Δt - 0.5
    else
        f = 0.0
        for j = 0:r.n
            f += binomial(r.n + j, j)*(r.n - j + 1)*(1 - Δt)^j
        end
        f*Δt^(r.n + 2)/(2r.n + 2)
    end
end

end
