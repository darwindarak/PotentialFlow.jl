module Freestreams

export Freestream

using ..Elements
using ..Motions

struct Freestream <: Element
    U::ComplexF64
end

Elements.@kind Freestream Singleton

Motions.induce_velocity(::ComplexF64, f::Freestream, t) = f.U
Motions.induce_velocity(f::Freestream, src, t) = nothing
Motions.induce_velocity!(vel, f::Freestream, src, t) = nothing

Elements.streamfunction(z::ComplexF64, f::Freestream) = -imag(f.U*z')
Elements.complexpotential(z::ComplexF64, f::Freestream) = f.U'*z

Elements.circulation(::Freestream) = 0.0
Elements.flux(::Freestream) = 0.0

Motions.self_induce_velocity!(vel, f::Freestream, t) = nothing

Motions.allocate_velocity(::Freestream) = ComplexF64

end
