import Base: *, \
struct ChebyshevTransform{T,I}
    dct!::FFTW.r2rFFTWPlan{T,(3,),true,1}
end

function plan_chebyshev_transform!(x::Vector{T}) where T
    ChebyshevTransform{T,true}(FFTW.plan_r2r!(x, FFTW.REDFT00))
end

function plan_chebyshev_transform(x::Vector{T}) where T
    ChebyshevTransform{T,false}(FFTW.plan_r2r!(x, FFTW.REDFT00))
end

(C::ChebyshevTransform{T,true}  * x::Vector{T}) where T = chebyshev_transform!(x, C.dct!)
(C::ChebyshevTransform{T,false} * x::Vector{T}) where T = chebyshev_transform(x, C.dct!)
(C::ChebyshevTransform{T,true}  \ A::Vector{T}) where T = inv_chebyshev_transform!(A, C.dct!)
(C::ChebyshevTransform{T,false} \ A::Vector{T}) where T = inv_chebyshev_transform(A, C.dct!)

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
    chebyshev_transform(x[, plan!])

Compute the discrete Chebyshev transform of `x`.

# Arguments

- `x`: samples of a function distributed along [extrema Chebyshev nodes](@ref chebyshev_nodes)
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function chebyshev_transform(x, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00))
    C = copy(x)
    chebyshev_transform!(C, plan!)
end

"""
    chebyshev_transform!(A,[ x, plan!])

Compute the discrete Chebyshev transform of `x` in-place.

# Arguments

- `A`: the output vector (also the input vector if `x` is omitted)
- `x`: an optional input vector with samples of a function distributed
  along [extrema Chebyshev nodes](@ref chebyshev_nodes)
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function chebyshev_transform!(x, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00))
    N = length(x)
    if N != length(plan!)
        throw(BoundsError("`x` must have the same size as the preplanned DCT"))
    end

    plan!*x

    s = 1/(N-1)

    @inbounds begin
        for n in 2:2:N
            x[n-1] *= s
            x[n]   *= -s
        end
        if isodd(N)
            x[N] *= s
        end

        x[1] /= 2
        x[N] /= 2
    end

    x
end

function chebyshev_transform!{T}(A::T, x::T, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00))
    copy!(A, x)
    chebyshev_transform!(A, plan!)
end

"""
    inv_chebyshev_transform(A[, plan!])

Perform the inverse Chebyshev transform.

This is the inverse of [`chebyshev_transform`](@ref).

# Fields

- `A`: the coefficients of a Chebyshev series
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function inv_chebyshev_transform(A, plan! = FFTW.plan_r2r!(A, FFTW.REDFT00))
    x = copy(A)
    inv_chebyshev_transform!(x, plan!)
end

"""
    inv_chebyshev_transform!(x, [C, plan!)

Perform the inverse Chebyshev transform in place.

This is the inverse of [`chebyshev_transform!`](@ref)

# Arguments

- `x`: the output vector to store the reconstructed function with
  points distributed along extrema Chebyshev nodes (also the input
  vector if `A` is omitted)
- `A`: an optional input vector with coefficients of a Chebyshev series
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function inv_chebyshev_transform!(A, plan! = FFTW.plan_r2r!(A, FFTW.REDFT00))
    N = length(A)
    if N != length(plan!)
        throw(BoundsError("`A` must have the same size as the preplanned DCT"))
    end

    @inbounds begin
        for n in 2:2:N-2
            A[n]   *= -0.5
            A[n+1] *=  0.5
        end
        if iseven(N)
            A[N] *= -1
        else
            A[N-1] *= -0.5
        end
    end

    plan!*A
end

function inv_chebyshev_transform!{T}(x::T, A::T, plan! = FFTW.plan_r2r!(A, FFTW.REDFT00))
    copy!(x, A)
    inv_chebyshev_transform!(x, plan!)
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
