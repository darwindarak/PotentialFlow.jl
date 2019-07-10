module Doublets

export Doublet

using ..Elements
using ..Motions

struct Doublet{T <: Number} <: Element
    z::ComplexF64
    S::T
end

Elements.kind(::Doublet) = Singleton
Elements.kind(::Type{Doublet{T}}) where T = Singleton

Elements.position(d::Doublet) = d.z
Elements.streamfunction(z::ComplexF64, d::Doublet) = imag(d.S/(z - d.z)/π)
Elements.complexpotential(z::ComplexF64, d::Doublet) = d.S/(z - d.z)/π


function Motions.induce_velocity(z::ComplexF64, d::Doublet, t)
    conj(-d.S/(z - d.z)^2/π)
end

end
