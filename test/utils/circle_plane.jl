# The functions in this file are meant to be used for testing the vortex
# evolution routines written using Chebyshev series.  They are not designed to
# fit with other rest of code.

struct JoukowskyMap
    c::ComplexF64
    α::Float64
end

function (m::JoukowskyMap)(ζ, derivative=0)
    if derivative == 0
        return joukowsky(ζ, m.c, m.α)
    elseif derivative == 1
        return joukowsky_d(ζ, m.α)
    elseif derivative == 2
        return joukowsky_d2(ζ, m.α)
    else
        throw(ArgumentError("JoukowskyMap is only implemented up to 2nd derivative"))
    end
end

joukowsky(ζ, c, α) = c + 0.5exp(im*α)*(ζ + 1/ζ)
joukowsky_d(ζ, α) = 0.5exp(im*α)*(1 - 1/ζ^2)
joukowsky_d2(ζ, α) = exp(im*α)/ζ^3

W_motion(ζ, ċ, α, Ω) = im*(imag(ċ*exp(-im*α))/ζ^2 + 0.5Ω/ζ^3)

function W_vortex(ζ, ζs, Γs)
    wᵥ = zero(Complex)
    for (ζᵥ, Γᵥ) in zip(ζs, Γs)
        wᵥ += Γᵥ*(1/(ζ - ζᵥ) - 1/(ζ - 1.0/conj(ζᵥ)))
    end
    return wᵥ/(2π*im)
end
