"""
    bound_circulation(plate[, s])

Compute the bound circulation on the plate.

`s` can be either a single normalized arc length coordinate (between
-1 and 1), or a whole array of coordinates.
"""
bound_circulation(plate) = bound_circulation(plate, plate.ss)

function bound_circulation(plate, s::Real)
    @get plate (A, α, B₀, B₁, L, Γ, N)

    γ = 0.0
    
    U₋ = 0.0
    U  = 1.0
    for n in 1:N-1
        γ += 2A[n]*U
        U₋, U = U, 2s*U - U₋
    end
    γ *= (s^2 - 1)

    γ += A[1] + 2Γ/(L*π)
    γ +=    2(A[0] - B₀)*s
    γ +=             -B₁*(2s^2 - 1)

    if abs2(γ) < 256eps()
        γ = 0.0
    else
        γ /= √(1 - s^2)
    end
    return γ
end

function bound_circulation(plate, ss::AbstractArray{T}) where {T <: Real}
    γs = zeros(Float64, length(ss))
    bound_circulation!(γs, plate, ss)
    return γs
end

"""
    bound_circulation!(γs, plate[, ss])

Compute the bound circulation of `plate` and store it in `γs`.

If an array, `ss`, with normalized arc length coordinates is omitted,
then the circulation will be computed at the plate's Chebyshev nodes.
"""
bound_circulation!(γs, plate) = bound_circulation!(γs, plate, plate.ss)

function bound_circulation!(γs, plate, ss)
    for (i, s) in enumerate(ss)
        γs[i] = bound_circulation(plate, s)
    end
    nothing
end

