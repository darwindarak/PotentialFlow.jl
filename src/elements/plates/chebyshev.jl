"""
    chebyshev_transform!{T}(C::T, x::T, plan! = FFTW.plan_r2r!(C, FFTW.REDFT00))

Compute the discrete Chebyshev transform of `x` and store the result in `C` (s = -1 → x[1], s = 1 → x[end]).
Also takes an optional `plan!` argument for a preplanned discrete cosine transform.
"""
function chebyshev_transform!{T}(C::T, x::T, plan! = FFTW.plan_r2r!(C, FFTW.REDFT00))
    N = length(x)
    for n in 1:N
        C[n] = x[n]/(N-1)
    end
    plan!*C

    for n in 2:2:N
        C[n] *= -1
    end

    C[1] /= 2
    C[N] /= 2
    nothing
end

"""
    chebyshev_transform(x)

Compute the discrete Chebyshev transform of `x` (s = -1 → x[1], s = 1 → x[end])
"""
function chebyshev_transform(x)
    C = zeros(eltype(x), length(x))
    chebyshev_transform!(C, x)
    return C
end

chebyshev_transform!(x, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00)) = chebyshev_transform!(x, x, plan!)

"""
    clenshaw_curtis_weights(N)

Returns the `N` weights used in the Clenshaw Curtis quadrature
"""
function clenshaw_curtis_weights(N)
    n = N - 1
    w = zeros(N)

    # Unlike the documentation, we move the 1/(N-1)
    # factor inside the r vector.
    w[1] = 1.0/n
    for k in 2:2:n
        w[k+1] = (1.0 + (-1)^k)/(1 - k^2)/n
    end

    # Scale by S
    w[1] *= 2; w[N] *= 2
    # Multiply by D
    FFTW.r2r!(w, FFTW.REDFT00)
    # Scale by S⁻¹
    w[1] /= 2; w[N] /= 2

    return w
end
