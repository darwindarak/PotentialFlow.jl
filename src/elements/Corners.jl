module Corners

export Corner, Wedge

using ..Elements
using ..Motions

"""
    Corner{N} <: Elements.Element

An immutable structure representing a corner (`N=1`)/wedge (`N=2`)

## Fields

- `σ::Real`: signed intensity
- `ν::Real`: interior angle factor (angle is νπ)
- `n₀::ComplexF64`: vector of bisector of interior angle
- `L::Real` : characteristic length of corner

## Constructors

- `Corner(σ,ν,θ)` : set up a corner with a bisector at angle `θ` and length `L=1`
- `Wedge(σ,ν,θ)`  : set up a wedge with a bisector at angle `θ` and length `L=1`
"""
struct Corner{N} <: Element
    σ::Real
    ν::Real
    n₀::ComplexF64
    L::Real
end

Corner(σ,ν,n₀,L) = Corner{1}(σ,ν,n₀,L)
Wedge(σ,ν,n₀,L) =  Corner{2}(σ,ν,n₀,L)

Corner(σ::Real,ν::Real,θ̄::Real) = Corner(σ,ν,exp(im*θ̄),1.0)
Wedge(σ::Real,ν::Real,θ̄::Real) = Wedge(σ,ν,exp(im*θ̄),1.0)


Elements.kind(::Corner) = Singleton
Elements.kind(::Type{Corner}) = Singleton

Base.angle(c::Corner) = angle(c.n₀)

function Elements.complexpotential(z::ComplexF64, c::Corner{N}) where {N}
    fact = exp(im*N*π/2)
    return fact*c.ν*c.σ*c.L^(1-1/c.ν).*(z*c.n₀').^(1/c.ν)
end

Elements.streamfunction(z::ComplexF64, c::Corner) = imag(complexpotential(z,c))

function Motions.induce_velocity(z::ComplexF64, c::Corner{N}, t) where {N}
    fact = exp(im*N*π/2)
    return conj(fact*c.σ*(z*c.n₀'/c.L).^(1/c.ν-1)*c.n₀')
end

end
