using Interpolations

function redistribute_points!(sheet, zs, Γs)
    sheet.Γs = Γs
    sheet.blobs = Vortex.Blob.(zs, compute_trapezoidal_weights(Γs), sheet.δ)
    sheet.zs = MappedPositions(Vortex.position, sheet.blobs, 0)
    sheet
end

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
    @assert 2 ≤ n ≤ length(sheet)-1
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
function remesh(sheet::Sheet, params::Tuple, Δs::Float64)
    L = arclengths(sheet)

    if L[end] < Δs
        warn("Cannot remesh, sheet length smaller then nominal remesh spacing")
        return Vortex.position.(sheet.blobs), sheet.Γs, L[end], params
    end

    knots = (L,)
    mode = Gridded(Linear())

    zspline = interpolate(knots, sheet.zs, mode)
    Γspline = interpolate(knots, sheet.Γs, mode)

    psplines = map(params) do p
        interpolate(knots, p, mode)
    end

    N = ceil(Int, L[end]/Δs)

    L₌ = linspace(0, L[end], N)

    z₌ = zspline[L₌]
    Γ₌ = Γspline[L₌]

    p₌ = map(psplines) do spline
        spline[L₌]
    end

    z₌, Γ₌, L[end], p₌

end
remesh(sheet::Sheet, Δs::Float64) = remesh(sheet, (), Δs)

function remesh!(sheet::Sheet, params::Tuple, Δs::Float64)
    z₌, Γ₌, L, p₌ = remesh(sheet, params, Δs)
    redistribute_points!(sheet, z₌, Γ₌)

    sheet, L, p₌
end
remesh!(sheet::Sheet, Δs::Float64) = remesh!(sheet, (), Δs)

function arclength(zs::AbstractVector)
    L = 0.0
    for i in 2:length(zs)
        L += abs(zs[i] - zs[i-1])
    end
    L
end
arclength(sheet::Sheet) = arclength(sheet.zs)

function arclengths(zs::AbstractVector)
    N = length(zs)

    Δs = zeros(N)
    for i in 2:N
        Δs[i] = abs(zs[i] - zs[i-1])
    end
    L = accumulate(+, Δs)
end
arclengths(sheet::Sheet) = arclengths(sheet.zs)

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
    Vortex.Sheets.filter!(sheet, Δf)

Apply Fourier filtering to the sheet position and strengths.  The
control points are redistributed to maintain a nominal point spacing
of of `Δs`, and the filtering removes any length scales smaller than
`Δf`.
"""
function filter!(sheet, Δs, Δf)
    z₌, Γ₌, l = remesh(sheet, Δs)
    if l > Δs
        filter_position!(z₌, Δf, l)
        redistribute_points!(sheet, z₌, Γ₌)
    else
        warn("Filter not applied, total sheet length smaller than nominal spacing")
        sheet
    end
end

function filter_position!(z₌::AbstractVector, Δf, L = arclength(z₌))
    Δs = abs(z₌[2] - z₌[1])

    @assert Δs < Δf

    # The in-place dct should be more performant, but it seems
    # to be giving different results on different trials.
    # Maybe this has to do with how FFTW plans the transforms?
    # dct!(z₌)
    ẑ = dct(z₌)

    cutoff = ceil(Int, 2L/Δf) + 1

    #z₌[cutoff:end] = zero(Complex128)
    ẑ[cutoff:end] = zero(Complex128)

    #idct!(z₌)
    idct(ẑ)
end

function filter_position!(sheet::Sheet, Δf, L = arclength(sheet))
    zs = Vortex.position.(sheet.blobs)
    filter_position!(z₌, Δf, L)

    redistribute_points!(sheet, z₌, sheet.Γs)
end
