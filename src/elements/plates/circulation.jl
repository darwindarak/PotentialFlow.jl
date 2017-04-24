function bound_circulation(s::Float64, plate)
    @get plate (A, ċ, α, α̇, L, Γ, N)

    B₀ = normal(ċ, α)
    B₁ = 0.5α̇*L

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

function bound_circulation!(γs, ss, plate)
    for (i, s) in enumerate(ss)
        γs[i] = bound_circulation(s, plate)
    end
    nothing
end

"""
    bound_circulation!(γs, plate)

Compute the bound circulation of the plate and store it in `γs`
"""
bound_circulation!(γs, plate) = bound_circulation!(γs, plate.ss, plate)

function bound_circulation(ss, plate)
    γs = zeros(Float64, length(ss))
    bound_circulation!(γs, ss, plate)
    return γs
end
bound_circulation(plate) = bound_circulation(plate.ss, plate)
