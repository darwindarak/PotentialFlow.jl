"""
    chebyshev_nodes(N, T = Float64)

Gives `N` extrema Chebyshev points of type `T` in [-1 ,1].
"""
function chebyshev_nodes(N, T = Float64)
    T[cos(θ) for θ in linspace(π, 0, N)]
end

"""
    eval_cheb(C, s)

Evaluate a Chebyshev series with coefficients `C` at point(s) `s`.

`s` can be either a single number or an array of numbers.
"""
eval_cheb(C, ss) = [eval_cheb(C, s) for s in ss]
function eval_cheb(C::AbstractArray{U}, s::S)::U where {U,S <: Real}
    -1 ≤ s ≤ 1 || throw(DomainError("s ∉ [-1,1]"))
    T₋ = s
    T  = one(S)

    f = zero(U)
    for c in C
        f += c*T
        T₋, T = T, 2s*T - T₋
    end
    f
end

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
    inv_chebyshev_transform!{T}(x::T, C::T, plan! = FFTW.plan_r2r!(C, FFTW.REDFT00))

Evaluate the Chebyshev series with coefficients `C` on the [extrema
nodes](@ref chebyshev_nodes) and store the results in `x`.

It is the inverse of [`chebyshev_transform!`](@ref), and similarly
takes an optional `plan!` argument for a preplanned discrete cosine
transform.
"""
function inv_chebyshev_transform!{T}(x::T, C::T, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00))
    N = length(C)
    x[1] = C[1]
    for n in 2:2:N-2
        x[n]   = -0.5C[n]
        x[n+1] =  0.5C[n+1]
    end
    if iseven(N)
        x[N] = -C[N]
    else
        x[N-1] = -0.5C[N-1]
        x[N] = C[N]
    end

    plan!*x
end
inv_chebyshev_transform!(C, plan! = FFTW.plan_r2r!(C, FFTW.REDFT00)) = inv_chebyshev_transform!(C, C, plan!)

"""
    inv_chebyshev_transform(C)

Inverts the coefficients `C` from a [discrete Chebyshev transform](@ref chebyshev_transform).
"""
function inv_chebyshev_transform(C)
    x = zeros(eltype(C), length(C))
    inv_chebyshev_transform!(x, C)
    return x
end

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
