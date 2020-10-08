module Sheets

export Sheet

using MappedArrays

using ..Blobs
using ..Elements
using Future: copy!
import ..Motions: position, mutually_induce_velocity!, self_induce_velocity!, advect!

const MappedPositions{T,R} = MappedArrays.ReadonlyMappedArray{ComplexF64,1,Array{Blob{T,R},1},typeof(Elements.position)} where {T,R}

mutable struct Sheet{T,R} <: Element
    blobs::Vector{Blob{T,R}}
    Ss::Vector{T}
    δ::Float64
    zs::MappedPositions{T,R}
end

Elements.unwrap(s::Sheet) = s.blobs

function Sheet(zs::AbstractVector{Complex{T}}, Ss::AbstractVector{S}, δ::Float64) where {T <: Number, S <: Number}
    dSs = compute_trapezoidal_weights(Ss)
    blobs = Blob{S,T}.(zs, dSs, δ)

    zs = mappedarray(Elements.position, blobs)

    Sheet{S,T}(blobs, copy(Ss), δ, zs)
end

function Sheet(blobs::Vector{Blob{T,R}}, Ss::AbstractVector, δ::Float64) where {T,R}
    newblobs = copy(blobs)
    zs = mappedarray(Elements.position, blobs)
    Sheet{T,R}(newblobs, copy(Ss), δ, zs)
end

Base.length(s::Sheet) = length(s.blobs)

function mutually_induce_velocity!(ws₁, ws₂, sheet₁::Sheet, sheet₂::Sheet, t)
    mutually_induce_velocity!(ws₁, ws₂, sheet₁.blobs, sheet₂.blobs, t)
    nothing
end

function self_induce_velocity!(ws, sheet::Sheet, t)
    self_induce_velocity!(ws, sheet.blobs, t)
end

function compute_trapezoidal_weights(Ss)
    N = length(Ss)

    dSs = similar(Ss)
    dSs[1] = 0.5*(Ss[2] - Ss[1])
    for i in 2:N-1
        dSs[i] = 0.5*(Ss[i+1] - Ss[i-1])
    end
    dSs[N] = 0.5*(Ss[N] - Ss[N-1])

    return dSs
end

function advect!(sheet₊::Sheet{S,R}, sheet₋::Sheet{S,R}, ws, Δt) where {S,R}
    if sheet₊ != sheet₋
        N₊ = length(sheet₊)
        N₋ = length(sheet₋)
        if N₊ != N₋
            resize!(sheet₊.blobs, N₋)
            resize!(sheet₊.Ss, N₋)
        end
        copy!(sheet₊.Ss, sheet₋.Ss)
        sheet₊.δ = sheet₋.δ
    end

    advect!(sheet₊.blobs, sheet₋.blobs, ws, Δt)
    sheet₊
end

include("sheets/surgery.jl")

end
