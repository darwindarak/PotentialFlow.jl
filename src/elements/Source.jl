module Source

import ..Points
import ..Blobs
import ..Elements: circulation, flux

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
const Point = Points.Point{ComplexF64}
(p::Point)(; z = p.z, S = imag(p.S)) = Point(z, S)

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
const Blob = Blobs.Blob{ComplexF64}
(b::Blob)(; z = b.z, S = imag(b.S), δ = b.δ) = Blob(z, S, δ)

function Base.show(io::IO, s::Blob)
    if iszero(real(s.S))
        print(io, "Source.Blob($(s.z), $(imag(s.S)), $(s.δ))")
    else
        print(io, "Blobs.Blob($(s.z), $(imag(s.S)), $(s.δ))")
    end
end
circulation(::Blob) = 0.0
flux(b::Blob) = imag(b.S)

end
