module Vortex

import ..Elements
import ..Elements: circulation, impulse

import ..Utils

import ..Points
import ..Blobs
import ..Sheets

#== Wrapper for a point vortex ==#

const Point = Points.Point{Float64}
Point(s::Point; z = s.z, S = s.S) = Point(z, S)

circulation(p::Point) = p.S
impulse(p::Point) = -im*p.z*p.S

Base.show(io::IO, s::Point) = print(io, "Vortex.Point($(s.z), $(s.S)")

#== Wrapper for a vortex blob ==#

const Blob = Blobs.Blob{Float64}
Blob(s::Blob; z = s.z, S = s.S, δ = s.δ) = Blob(z, S, δ)

circulation(b::Blob) = b.S
impulse(b::Blob) = -im*b.z*b.S
Base.show(io::IO, s::Blob) = print(io, "Vortex.Blob($(s.z), $(s.S), $(s.δ)")

#== Wrapper for a vortex sheet ==#

const Sheet = Sheets.Sheet{Float64}
Sheet(blobs::Vector{Blob}, Ss::AbstractVector{Float64}, δ::Float64) = Sheets.Sheet(blobs, Ss, δ)
Sheet(zs::AbstractVector,  Ss::AbstractVector{Float64}, δ::Float64) = Sheets.Sheet(zs, Ss, δ)
function Base.show(io::IO, s::Sheet)
    L = Sheets.arclength(s)
    print(io, "Vortex Sheet: L ≈ $(round(L, 3)), Γ = $(round(s.Ss[end] - s.Ss[1], 3)), δ = $(round(s.δ, 3))")
end

circulation(s::Sheet) = s.Ss[end] - s.Ss[1]
impulse(s::Sheet) = impulse(s.blobs)

end
