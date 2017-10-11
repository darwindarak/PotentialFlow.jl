"""
    strength(plate[, s])

Compute the bound vortex sheet strength of the plate.

`s` can be either a single normalized arc length coordinate (between
-1 and 1), or a whole array of coordinates.
"""
strength(plate) = strength(plate, plate.ss)

function strength(plate, s::Real)
    @get plate (A, α, B₀, B₁, L, Γ, N)

    γ = 0.0

    U₋ = 1.0
    U  = 2s
    for n in 2:N-1
        γ += 2A[n]*U
        U₋, U = U, 2s*U - U₋
    end
    γ *= (s^2 - 1)

    γ += 2Γ/(L*π)
    γ += 2(A[0] - B₀)*s
    γ +=  (A[1] - B₁)*(2s^2 - 1)

    if abs2(γ) < 256eps()
        γ = 0.0
    else
        γ /= √(1 - s^2)
    end
    return γ
end

function strength(plate, ss::AbstractArray{T}) where {T <: Real}
    γs = zeros(Float64, length(ss))
    strength!(γs, plate, ss)
    return γs
end

"""
    strength!(γs, plate[, ss])

Compute the bound vortex sheet strength of `plate` and store it in `γs`.

If an array, `ss`, with normalized arc length coordinates is omitted,
then the strength will be computed at the plate's Chebyshev nodes.
"""
function strength!(γs, plate, ss = plate.ss)
    for (i, s) in enumerate(ss)
        γs[i] = strength(plate, s)
    end
    nothing
end

"""
    bound_circulation(plate[, s])

Compute the bound circulation between the trailing edge of the plate to `s`.

`s` can be either a single normalized arc length coordinate (between
-1 and 1), or a whole array of coordinates.
"""
bound_circulation(plate) = bound_circulation(plate, plate.ss)

function bound_circulation(plate, s::Real)
    @get plate (A, B₀, B₁, L, Γ)
    _bound_circulation(A, B₀, B₁, L, Γ, s)
end

function _bound_circulation(A, B₀, B₁, L, Γ, s)
    N = length(A)

    Γₛ = 0.0


    Γₛ  = 2(A[0] - B₀)
    Γₛ +=  (A[1] - B₁)*s

    U₋₂ = 1.0
    U₋  = 2s
    U   = 4s^2 - 1

    for n in 2:N-1
        Γₛ += A[n]*(U/(n+1) - U₋₂/(n-1))
        U₋₂, U₋, U = U₋, U, 2s*U - U₋
    end
    Γₛ *= -0.5L*√(1 - s^2)

    Γₛ -= Γ*(acos(s)/π - 1)

    return Γₛ
end

function bound_circulation(plate, ss::AbstractArray{T}) where {T <: Real}
    Γs = zeros(Float64, length(ss))
    bound_circulation!(Γs, plate, ss)
    return Γs
end

"""
    bound_circulation!(Γs, plate[, ss])

Compute the bound circulation between the trailing edge of the plate to `ss`, then store it in `Γs`.

If an array, `ss`, with normalized arc length coordinates is omitted,
then the circulation will be computed at the plate's Chebyshev nodes.
"""
function bound_circulation!(Γs, plate, ss = plate.ss)
    for (i, s) in enumerate(ss)
        Γs[i] = bound_circulation(plate, s)
    end
    nothing
end
