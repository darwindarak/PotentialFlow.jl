using Interpolations
using FFTW

"""
    Sheets.redistribute_points!(sheet, zs, Γs)

Returns the modified sheet with replacement control points at positions `zs` and strength `Γs`.

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> sys = (sheet,)
(Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2,)

julia> Sheets.redistribute_points!(sheet, 0:0.2:2, 0.0:0.5:5)
Vortex Sheet: L ≈ 2.0, Γ = 5.0, δ = 0.2

julia> sys[1]
Vortex Sheet: L ≈ 2.0, Γ = 5.0, δ = 0.2
```
"""
function redistribute_points!(sheet::Sheet{T}, zs, Ss) where T
    if !(sheet.Ss === Ss)
        resize!(sheet.Ss, length(Ss))
        copy!(sheet.Ss, Ss)
    end

    sheet.blobs = Blob{T}.(zs, compute_trapezoidal_weights(Ss), sheet.δ)
    sheet.zs = mappedarray(Elements.position, sheet.blobs)
    sheet
end

"""
    Sheets.split!(sheet, n::Int)

Remove segments `0:n` from `sheet`, and return those segments as a new sheet.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> sheet₋ = Sheets.split!(sheet, 5)
Vortex Sheet: L ≈ 0.4, Γ = 4.0, δ = 0.2

julia> sheet
Vortex Sheet: L ≈ 0.6, Γ = 6.0, δ = 0.2
```
"""
function split!(sheet::Sheet{T}, n::Int) where T
    @assert 2 < n < length(sheet) - 2

    blobs = splice!(sheet.blobs, 1:n-1)
    Ss = splice!(sheet.Ss, 1:n-1)
    push!(Ss, sheet.Ss[1])
    push!(blobs, Blob{T}(sheet.blobs[1].z, 0.5(Ss[end] - Ss[end-1]), sheet.δ))

    sheet.blobs[1] = Blob{T}(sheet.blobs[1].z, 0.5(sheet.Ss[2] - sheet.Ss[1]), sheet.δ)

    Sheet(blobs, Ss, sheet.δ)
end

"""
    Sheets.truncate!(sheet, n::Int)

Remove segments `0:n` from `sheet`, and return the circulation in those segments.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> Sheets.truncate!(sheet, 5)
4.0
```
"""
function truncate!(sheet::Sheet{T}, n::Int) where T
    @assert 2 ≤ n ≤ length(sheet)-1
    ΔS = sheet.Ss[n] - sheet.Ss[1]
    splice!(sheet.Ss, 1:n-1)
    splice!(sheet.blobs, 1:n-1)
    sheet.blobs[1] = Blob{T}(sheet.blobs[1].z, 0.5(sheet.Ss[2] - sheet.Ss[1]), sheet.δ)
    return ΔS
end

"""
    Sheets.remesh(sheet, Δs::Float64 , params::Tuple = ())

Uniformly redistribute the control points of the sheet to have a nominal spacing of `Δs`.
Material quantities that should be redistributed along with the control points can be passed in as elements of `params`.

Returns the tuple `(z₌, Γ₌, L [, p₌])` where

- `z₌` is an array with the positions of the uniformly distributed points
- `Γ₌` is circulation interpolated onto `z₌`
- `L` is total length of the sheet
- `p₌` is a tuple containing the material quantities from `params` interpolated onto `z₌`

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> age = collect(10.0:-1:0);

julia> Sheets.remesh(sheet, 0.2, (age, ))
(Complex{Float64}[0.0+0.0im, 0.25+0.0im, 0.5+0.0im, 0.75+0.0im, 1.0+0.0im], [0.0, 2.5, 5.0, 7.5, 10.0], 1.0, ([10.0, 7.5, 5.0, 2.5, 0.0],))
```
"""
function remesh(sheet::Sheet{S}, Δs::Float64, params::Tuple = ()) where S
    L = arclengths(sheet)

    if L[end] < Δs
        @warn("Cannot remesh, sheet length smaller than nominal spacing")
        return Elements.position.(sheet.blobs), sheet.Ss, L[end], params
    end

    knots = (L,)
    mode = Gridded(Linear())

    zspline = interpolate(knots, sheet.zs, mode)
    Γspline = interpolate(knots, sheet.Ss, mode)

    psplines = map(params) do p
        interpolate(knots, p, mode)
    end

    N = ceil(Int, L[end]/Δs)

    L₌ = range(0, L[end], length=N)

    z₌ = zspline(L₌)
    S₌ = Γspline(L₌)

    p₌ = map(psplines) do spline
        spline(L₌)
    end

    z₌, S₌, L[end], p₌
end

"""
    Sheets.remesh!(sheet::Sheet, Δs::Float64, params::Tuple = ())

Same as [`Sheets.remesh`](@ref), except `sheet` is replaced
internally by a uniformly interpolated control points.
Returns the tuple (sheet, L, p₌) where

- `sheet` is the modified sheet
- `L` is total length of the sheet
- `p₌` is a tuple containing the material quantities from `params` interpolated onto the new control points of `sheet`

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> age = collect(10.0:-1:0);

julia> Sheets.remesh!(sheet, 0.2, (age,));

julia> Elements.position.(sheet.blobs)
5-element Array{Complex{Float64},1}:
  0.0 + 0.0im
 0.25 + 0.0im
  0.5 + 0.0im
 0.75 + 0.0im
  1.0 + 0.0im

julia> age
5-element Array{Float64,1}:
 10.0
  7.5
  5.0
  2.5
  0.0
```
"""
function remesh!(sheet::Sheet{S}, Δs::Float64, params::Tuple = ()) where S
    z₌, S₌, L, p₌ = remesh(sheet, Δs, params)
    if L < Δs
        return sheet, L, params
    end
    redistribute_points!(sheet, z₌, S₌)

    for i in 1:length(params)
        resize!(params[i], length(p₌[i]))
        copy!(params[i], p₌[i])
    end

    sheet, L, params
end

"""
    Sheets.arclength(s)

Compute the polygonal arc length of `s`, where `s` can be either an
vector of complex numbers or a `Vortex.Sheet`.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> Sheets.arclength(sheet)
1.0
"""
function arclength(zs::AbstractVector)
    L = 0.0
    for i in 2:length(zs)
        L += abs(zs[i] - zs[i-1])
    end
    L
end
arclength(sheet::Sheet) = arclength(sheet.zs)

"""
    Sheets.arclengths(s)

Cumulative sum of the polygonal arc length of `s`, where `s` can be
either an vector of complex numbers or a `Vortex.Sheet`.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> Sheets.arclengths(sheet)
11-element Array{Float64,1}:
 0.0
 0.1
 0.2
 0.3
 0.4
 0.5
 0.6
 0.7
 0.8
 0.9
 1.0
```
"""
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
    Sheets.append_segment!(sheet::Sheet, z, Γ)

Append a new segment with circulation `Γ` extending from the end of the sheet to `z`.

# Example

```jldoctest
julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)
Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2

julia> sheet.blobs[end]
Vortex.Blob(1.0 + 0.0im, 0.5, 0.2)

julia> Sheets.append_segment!(sheet, 1.1, 2.0)

julia> sheet
Vortex Sheet: L ≈ 1.1, Γ = 12.0, δ = 0.2

julia> sheet.blobs[end]
Vortex.Blob(1.1 + 0.0im, 1.0, 0.2)
```
"""
function append_segment!(sheet::Sheet{T}, z, S::T) where T
    b₋ = sheet.blobs[end]
    sheet.blobs[end] = Blob{T}(b₋.z, b₋.S + 0.5S, b₋.δ)
    push!(sheet.blobs, Blob{T}(z, 0.5S, sheet.δ))
    push!(sheet.Ss, sheet.Ss[end] + S)
    nothing
end

"""
    Sheets.filter!(sheet, Δs, Δf[, params])

Redistribute and filter the control points of a vortex sheet

# Arguments

- `sheet`: the vortex sheet to be modified
- `Δs`: the nominal spacing between the uniform points
- `Δf`: the minimum length scale that the filter should allow to pass
  through
- `params`: an optional tuple of vectors containing material properties

# Returns

If `params` is passed in, then its vectors will be overwritten by
their interpolated values on the new control points, and the function
returns the tuple (sheet, params).
Otherwise, it returns (sheet, ())
"""
function filter!(sheet, Δs, Δf, params::Tuple = ())
    z₌, S₌, L, p₌ = remesh(sheet, Δs, params)

    for i in 1:length(params)
        resize!(params[i], length(p₌[i]))
        copy!(params[i], p₌[i])
    end

    if L > Δs
        filter_position!(z₌, Δf, L)
        redistribute_points!(sheet, z₌, S₌), params
    else
        @warn("Filter not applied, total sheet length smaller than nominal spacing")
        sheet, params
    end
end

"""
    Sheets.filter_position!(s, Δf, L = arclength(z₌))

Filter out any length scales in `s` that is smaller than `Δf`, storing the result back in `s`.
`s` can be either a vector of complex positions, or a `Vortex.Sheet`.
"""
function filter_position!(z₌::AbstractVector{T}, Δf, L = arclength(z₌)) where T
    Δs = abs(z₌[2] - z₌[1])

    @assert Δs < Δf

    F = plan_dct!(z₌)
    F * z₌

    cutoff = ceil(Int, 2L/Δf) + 1

    z₌[cutoff:end] .= zero(ComplexF64)

    # F \ z₌
    # Out-of-date syntax in FFTW for planned inverse transforms
    # manually construct inverse transform for now
    F⁻¹ = let X = Array{T}(undef, F.plan.sz),
             iK = FFTW.inv_kind[FFTW.REDFT10]
           FFTW.DCTPlan{T,iK,true}(FFTW.plan_r2r!(X, iK, F.region, flags=F.plan.flags),
                              F.r, F.nrm, F.region)
       end

    F⁻¹ * z₌
end

function filter_position!(sheet::Sheet, Δf, L = arclength(sheet))
    zs = Elements.position.(sheet.blobs)
    filter_position!(zs, Δf, L)

    redistribute_points!(sheet, zs, sheet.Ss)
end
