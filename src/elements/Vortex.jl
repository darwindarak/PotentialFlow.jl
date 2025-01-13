module Vortex

import ..Elements
import ..Elements: circulation, flux, impulse, angularimpulse, seed_zeros, seed_position,
                    seed_strength, kind, seed_position!, seed_strength!

import ..Utils: dualize, ComplexGradientConfig, seed!, Partials

import ..Points
import ..Blobs
import ..Sheets

#== Wrapper for a point vortex ==#

"""
    Vortex.Point(z::ComplexF64, Γ::Float64)

A point vortex located at `z` with circulation `Γ`.

A new point vortex can be created from an existing one by treating the
existing point vortex as a function and passing in the parameter you
want to change as keyword arguments.
For example,
```jldoctest
julia> p = Vortex.Point(1.0, 1.0)
Vortex.Point(1.0 + 0.0im, 1.0)

julia> p()
Vortex.Point(1.0 + 0.0im, 1.0)

julia> p(Γ = 2.0)
Vortex.Point(1.0 + 0.0im, 2.0)
```
""" Vortex.Point(::ComplexF64,::Float64)
#const Point = Points.Point{Float64,Float64}
const Point = Points.Point{T,R} where {T<: Real, R<:Real}
Point(z::Complex{R},Γ::T;period=Inf) where {T<:Real,R<:Real} = Points.Point{T}(z,Γ,period)
Point(z::Real,Γ::T;period=Inf) where {T} = Points.Point{T}(complex(z),Γ,period)

(p::Point)(; z = p.z, Γ = p.S, period=p.period) = Point(z, Γ;period=period)

circulation(p::Point) = p.S
flux(::Point) = 0.0
impulse(p::Point) = -im*p.z*p.S
angularimpulse(p::Point) = -0.5*p.z*conj(p.z)*p.S

Base.show(io::IO, s::Point) = print(io, "Vortex.Point($(s.z), $(s.S))")


#== Wrapper for a vortex blob ==#

"""
    Vortex.Blob(z::ComplexF64, Γ::Float64, δ::Float64)

A regularized point vortex located at `z` with circulation `Γ` and blob radius `δ`.

A new vortex blob can be created from an existing one by treating the
existing blob as a function and passing in the parameter you want to
change as keyword arguments.
For example,
```jldoctest
julia> b = Vortex.Blob(1.0, 1.0, 0.1)
Vortex.Blob(1.0 + 0.0im, 1.0, 0.1)

julia> b()
Vortex.Blob(1.0 + 0.0im, 1.0, 0.1)

julia> b(Γ = 2.0, δ = 0.01)
Vortex.Blob(1.0 + 0.0im, 2.0, 0.01)
```
""" Vortex.Blob(::ComplexF64,::Float64,::Float64)
#const Blob = Blobs.Blob{Float64,Float64}
const Blob = Blobs.Blob{T,R,P} where {T<:Real,R<:Real,P}
Blob(z::Complex{R},Γ::T,δ::Float64;period=Inf) where {T<:Real,R<:Real} = Blobs.Blob{T}(z,Γ,δ,period)
Blob(z::Real,Γ::T,δ;period=Inf) where {T<:Real} = Blobs.Blob{T}(complex(z),Γ,δ,period)



(b::Blob)(; z = b.z, Γ = b.S, δ = b.δ, period=b.period) = Blob(z, Γ, δ; period = period)

circulation(b::Blob) = b.S
flux(::Blob) = 0.0
impulse(b::Blob) = -im*b.z*b.S
angularimpulse(b::Blob) = -0.5*b.z*conj(b.z)*b.S
Base.show(io::IO, s::Blob) = print(io, "Vortex.Blob($(s.z), $(s.S), $(s.δ), $(s.period))")



#== Wrapper for a vortex sheet ==#

"""
    Vortex.Sheet <: Elements.Element

A vortex sheet represented by vortex blob control points

## Fields

- `blobs`: the underlying array of vortex blobs
- `Ss`: the cumulated sum of circulation starting from the first control point
- `δ`: the blob radius of all the vortex blobs
- `zs`: a mapped array that accesses the position of each control point

## Constructors:

- `Vortex.Sheet(blobs, Γs, δ)`
- `Vortex.Sheet(zs, Γs, δ)` where `zs` is an array of positions for the control points
"""
const Sheet = Sheets.Sheet{T,R,P} where {T<:Real,R<:Real,P}
Sheet(blobs::Vector{Blob{T,R,P}}, Ss::AbstractVector{Float64}, δ::Float64) where {T,R,P }= Sheets.Sheet(blobs, Ss, δ)
Sheet(zs::AbstractVector,  Ss::AbstractVector{Float64}, δ::Float64;period=Inf) = Sheets.Sheet(zs, Ss, δ)

function Base.show(io::IO, s::Sheet)
    L = Sheets.arclength(s)
    print(io, "Vortex Sheet: L ≈ $(round(L, digits=3)), Γ = $(round(s.Ss[end] - s.Ss[1], digits=3)), δ = $(round(s.δ, digits=3))")
end

circulation(s::Sheet) = s.Ss[end] - s.Ss[1]

include("vortex/diff.jl")

end
