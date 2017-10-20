module Doublets

using ..Elements
using ..Motions

struct Doublet{T <: Number} <: Element
    z::Complex128
    S::T
end

Elements.kind(::Doublet) = Singleton
Elements.kind(::Type{Doublet{T}}) where T = Singleton

Elements.position(d::Doublet) = d.z
Elements.streamfunction(z::Complex128, d::Doublet) = imag(d.S/(z - d.z)/π)

function Motions.induce_velocity(z::Complex128, d::Doublet, t)
    conj(-d.S/(z - d.z)^2)
end

end
