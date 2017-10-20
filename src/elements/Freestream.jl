struct Freestream <: Element
    U::Complex128
end

Elements.@kind Freestream Singleton

Motions.induce_velocity(::Complex128, f::Freestream, t) = f.U
Motions.induce_velocity(f::Freestream, src, t) = nothing
Motions.induce_velocity!(vel, f::Freestream, src, t) = nothing

Elements.streamfunction(z::Complex128, f::Freestream) = imag(f.U*z)

Motions.self_induce_velocity!(vel, f::Freestream, t) = nothing

Motions.allocate_velocity(::Freestream) = nothing
