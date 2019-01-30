module Chebyshev

using FFTW
using Future: copy!

import Base: *, \
struct Transform{T,I}
    dct!::FFTW.r2rFFTWPlan{T,(3,),true,1}
end

function plan_transform!(x::Vector{T}) where T
    Transform{T,true}(FFTW.plan_r2r!(x, FFTW.REDFT00))
end

function plan_transform(x::Vector{T}) where T
    Transform{T,false}(FFTW.plan_r2r!(x, FFTW.REDFT00))
end

(C::Transform{T,true}  * x::Vector{T}) where T = transform!(x, C.dct!)
(C::Transform{T,false} * x::Vector{T}) where T = transform(x, C.dct!)
(C::Transform{T,true}  \ A::Vector{T}) where T = inv_transform!(A, C.dct!)
(C::Transform{T,false} \ A::Vector{T}) where T = inv_transform(A, C.dct!)

"""
    Chebyshev.nodes(N, T = Float64)

Gives `N` extrema Chebyshev points of type `T` in [-1 ,1].
"""
function nodes(N, T = Float64)
    T[cos(θ) for θ in range(π, 0, length=N)]
end

"""
    Chebyshev.firstkind(C, s[, offset = 0])

Evaluate a Chebyshev series of the first kind.

# Arguments

- `C`: a vector with the coefficients of the series
- `s`: a scalar or an array of evaluation points ∈ [-1, 1]
- `offset`: an optional starting index of the series summation
"""
firstkind(A, ss, offset = 0) = [firstkind(A, s, offset) for s in ss]
function firstkind(A::AbstractArray{C}, s::S, offset = 0)::C where {C,S <: Real}
    -1 ≤ s ≤ 1 || throw(DomainError("s ∉ [-1,1]"))
    T₋ = s
    T  = one(S)

    for i in 1:offset
        T₋, T = T, 2s*T - T₋
    end

    N = length(A)
    f = zero(C)
    for i in 1+offset:N
        f += A[i]*T
        T₋, T = T, 2s*T - T₋
    end
    f
end

"""
    Chebyshev.secondkind(C, s[, offset = 0])

Evaluate a Chebyshev series of the second kind.

# Arguments

- `C`: a vector with the coefficients of the series
- `s`: a scalar or an array of evaluation points ∈ [-1, 1]
- `offset`: an optional starting index of the series summation
"""
secondkind(C, ss, offset = 0) = [secondkind(C, s, offset) for s in ss]
function secondkind(A::AbstractArray{C}, s::S, offset = 0)::C where {C,S <: Real}
    -1 ≤ s ≤ 1 || throw(DomainError("s ∉ [-1,1]"))
    U₋ = zero(s)
    U  = one(S)

    for i in 1:offset
        U₋, U = U, 2s*U - U₋
    end

    N = length(A)
    f = zero(C)
    for i in 1+offset:N
        f += A[i]*U
        U₋, U = U, 2s*U - U₋
    end
    f
end

"""
    Chebyshev.transform(x[, plan!])

Compute the discrete Chebyshev transform of `x`.

# Arguments

- `x`: samples of a function distributed along [extrema Chebyshev nodes](@ref Chebyshev.nodes)
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function transform(x, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00))
    C = copy(x)
    transform!(C, plan!)
end

"""
    Chebyshev.transform!(A,[ x, plan!])

Compute the discrete Chebyshev transform of `x` in-place.

# Arguments

- `A`: the output vector (also the input vector if `x` is omitted)
- `x`: an optional input vector with samples of a function distributed
  along [extrema Chebyshev nodes](@ref Chebyshev.nodes)
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function transform!(x, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00))
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

function transform!(A::T, x::T, plan! = FFTW.plan_r2r!(x, FFTW.REDFT00)) where T
    copy!(A, x)
    transform!(A, plan!)
end

"""
    Chebyshev.inv_transform(A[, plan!])

Perform the inverse Chebyshev transform.

This is the inverse of [`Chebyshev.transform`](@ref).

# Fields

- `A`: the coefficients of a Chebyshev series
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function inv_transform(A, plan! = FFTW.plan_r2r!(A, FFTW.REDFT00))
    x = copy(A)
    inv_transform!(x, plan!)
end

"""
    Chebyshev.inv_transform!(x, [C, plan!)

Perform the inverse Chebyshev transform in place.

This is the inverse of [`Chebyshev.transform!`](@ref)

# Arguments

- `x`: the output vector to store the reconstructed function with
  points distributed along extrema Chebyshev nodes (also the input
  vector if `A` is omitted)
- `A`: an optional input vector with coefficients of a Chebyshev series
- `plan!`: an optional pre-planned in-place DCT-I used to compute the transform
"""
function inv_transform!(A, plan! = FFTW.plan_r2r!(A, FFTW.REDFT00))
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

function inv_transform!(x::T, A::T, plan! = FFTW.plan_r2r!(A, FFTW.REDFT00)) where T
    copy!(x, A)
    Chebyshev.inv_transform!(x, plan!)
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

end
