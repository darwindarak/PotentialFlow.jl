using Interpolations

"""
    Vortex.Sheets.split!(sheet, n::Int)

Remove segments `0:n` from `sheet`, and return those segments as a new sheet.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> sheet₋ = Vortex.Sheets.split!(sheet, 5)
Vortex Sheet: L ≈ 0.4, Γ = 4.0, δ = 0.2

julia> sheet
Vortex Sheet: L ≈ 0.6, Γ = 6.0, δ = 0.2
```
"""
function split!(sheet, n::Int)
    @assert 2 < n < length(sheet) - 2

    blobs = splice!(sheet.blobs, 1:n-1)
    Γs = splice!(sheet.Γs, 1:n-1)
    push!(Γs, sheet.Γs[1])
    push!(blobs, Vortex.Blob(sheet.blobs[1].z, 0.5(Γs[end] - Γs[end-1]), sheet.δ))

    sheet.blobs[1] = Vortex.Blob(sheet.blobs[1].z, 0.5(sheet.Γs[2] - sheet.Γs[1]), sheet.δ)

    Vortex.Sheet(blobs, Γs, sheet.δ)
end

"""
    Vortex.Sheets.truncate!(sheet, n::Int)

Remove segments `0:n` from `sheet`, and return the circulation in those segments.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> Vortex.Sheets.truncate!(sheet, 5)
4.0
```
"""
function truncate!(sheet, n::Int)
    @assert 2 ≤ n < length(sheet)-1
    ΔΓ = sheet.Γs[n] - sheet.Γs[1]
    splice!(sheet.Γs, 1:n-1)
    splice!(sheet.blobs, 1:n-1)
    sheet.blobs[1] = Vortex.Blob(sheet.blobs[1].z, 0.5(sheet.Γs[2] - sheet.Γs[1]), sheet.δ)
    return ΔΓ
end

"""
    Vortex.Sheets.remesh!(sheet, zs, Γs)

Redistribute the control points of the sheet to lie on `zs` with circulation `Γs`.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> Vortex.Sheets.remesh!(sheet, 0:0.2:2, 2sheet.Γs)
Vortex Sheet: L ≈ 2.0, Γ = 20.0, δ = 0.2
```
"""
function remesh!(sheet, zs, Γs)
    sheet.Γs = Γs
    sheet.blobs = Vortex.Blob.(zs, compute_trapezoidal_weights(Γs), sheet.δ)
    sheet
end

"""
    Vortex.Sheets.append_segment!(sheet::Sheet, z, Γ)

Append a new segment with circulation `Γ` extending from the end of the sheet to `z`.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> Vortex.append_segment!(sheet, 1.1, 2.0)
Vortex Sheet: L ≈ 1.1, Γ = 12.0, δ = 0.2
```
"""
function append_segment!(sheet::Sheet, z, Γ)
    b₋ = sheet.blobs[end]
    sheet.blobs[end] = Vortex.Blob(b₋.z, b₋.Γ + 0.5Γ, b₋.δ)
    push!(sheet.blobs, Vortex.Blob(z, 0.5Γ, sheet.δ))
    push!(sheet.Γs, sheet.Γs[end] + Γ)
    nothing
end

"""
    Vortex.Sheets.filter!(sheet, Δs, Δf)

Apply Fourier filtering to the sheet position and strengths.  The
control points are redistributed to maintain a nominal point spacing
of of `Δs`, and the filtering removes any length scales smaller than
`Δf`.
"""
function filter!(sheet, Δs, Δf)
    zs = Vortex.position.(sheet.blobs)
    zs, (Γs,) = filter_by_arclength(zs, Δs, Δf, sheet.Γs)
    remesh!(sheet, zs, Γs)
end

function filter_by_arclength(z, Δs, Δf, properties...)
    @assert Δs < Δf

    l = vcat(0, accumulate(+, abs.(diff(z))))

    n = round(Int, l[end]/Δs)
    if n ≤ 1
        warn("Sampling interval too large for arc length")
        return z, properties
    end

    z_spline = interpolate((l,), z, Gridded(Linear()))

    splines = map(properties) do p
        interpolate((l,), p, Gridded(Linear()))
    end

    L = linspace(0, l[end], n)
    ΔL = step(L)

    z₌ = z_spline[L]
    splines₌ = map(splines) do spline
        spline[L]
    end

    ẑ = dct(z₌)
    Ŝ = dct.(splines₌)

    cutoff = ceil(Int, 2L[end]/Δf) + 1
    ẑ[cutoff:end] = zero(Complex128)
    foreach(Ŝ) do S
        S[cutoff:end] = 0.0
    end

    idct(ẑ), idct.(Ŝ)
end
