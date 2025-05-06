module Source

import ..Points
import ..Blobs
import ..Elements: circulation, flux, kind, seed_zeros, seed_position, seed_strength,
                    seed_position!, seed_strength!

import ..Utils: dualize, ComplexGradientConfig, seed!, Partials


#== Wrapper for a point source ==#

"""
    Source.Point(z::ComplexF64, S::Float64)

A point source located at `z` with strength `S`.

A new point source can be created from an existing one by treating the
existing source as a function and passing in the parameter you
want to change as keyword arguments.
For example,
```jldoctest
julia> p = Source.Point(1.0, 1.0)
Source.Point(1.0 + 0.0im, 1.0)

julia> p()
Source.Point(1.0 + 0.0im, 1.0)

julia> p(S = 2.0)
Source.Point(1.0 + 0.0im, 2.0)
```
"""
const Point = Points.Point{T,R} where {T<:Complex, R <: Real}
Point(z::Complex{R},S::T;period=Inf) where {T<:Real,R<:Real} = Points.Point{Complex{T}}(z,S,period)
Point(z::Real,S::T;period=Inf) where {T} = Points.Point{Complex{T}}(complex(z),S,period)

(p::Point)(; z = p.z, S = imag(p.S), period=p.period) = Point(z, S; period=period)

function Base.show(io::IO, s::Point)
    if iszero(real(s.S))
        print(io, "Source.Point($(s.z), $(imag(s.S)))")
    else
        print(io, "Points.Point($(s.z), $(s.S))")
    end
end
flux(p::Point) = imag(p.S)
circulation(::Point) = 0.0


#== Wrapper for a blob source ==#

"""
    Source.Blob(z::ComplexF64, S::Float64, δ::Float64)

A regularized point source located at `z` with strength `S` and blob radius `δ`.

A new blob source can be created from an existing one by treating the
existing blob as a function and passing in the parameter you want to
change as keyword arguments.
For example,
```jldoctest
julia> b = Source.Blob(1.0, 1.0, 0.1)
Source.Blob(1.0 + 0.0im, 1.0, 0.1)

julia> b()
Source.Blob(1.0 + 0.0im, 1.0, 0.1)

julia> b(S = 2.0, δ = 0.01)
Source.Blob(1.0 + 0.0im, 2.0, 0.01)
```
"""
const Blob = Blobs.Blob{T,R,P} where {T<:Complex, R <: Real, P}
Blob(z::Complex{R},S::T,δ::Float64;period=Inf) where {T<:Real,R<:Real} = Blobs.Blob{Complex{T}}(z,S,δ,period)
Blob(z::Real,S::T,δ;period=Inf) where {T<:Real} = Blobs.Blob{Complex{T}}(complex(z),S,δ,period)


(b::Blob)(; z = b.z, S = imag(b.S), δ = b.δ, period = b.period) = Blob(z, S, δ; period=period)

function Base.show(io::IO, s::Blob)
    if iszero(real(s.S))
        if isinf(s.period)
            print(io, "Source.Blob($(s.z), $(imag(s.S)), $(s.δ))")
        else                                            
            print(io, "Source.Blob($(s.z), $(imag(s.S)), $(s.δ), $(s.period))")
        end
    else
        if isinf(s.period)
            print(io, "Blobs.Blob($(s.z), $(imag(s.S)), $(s.δ))")
        else
            print(io, "Blobs.Blob($(s.z), $(imag(s.S)), $(s.δ), $(s.period))")
        end
    end
end
circulation(::Blob) = 0.0
flux(b::Blob) = imag(b.S)


include("source/diff.jl")

end
