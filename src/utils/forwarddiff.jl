
using ForwardDiff
using DiffRules

import ForwardDiff:value,partials,derivative,extract_derivative, Dual, seed!
const AMBIGUOUS_TYPES = (AbstractFloat, Irrational, Integer, Rational, Real, RoundingMode, ComplexF64)

# ComplexComplexDual is designated for derivatives of complex functions
# with respect to complex numbers
# ComplexRealDual is for derivatives of complex functions with respect to
# real numbers
# ComplexDual encompasses all complex dual types
const ComplexDual{T,V,N} = Complex{Dual{T,V,N}}
const ComplexComplexDual{T,V} = ComplexDual{T,V,2}
const ComplexRealDual{T,V} = ComplexDual{T,V,1}
const RealComplexDual{T,V} = Dual{T,V,2}


@inline function ComplexComplexDual{T}(z::Number,dz::Number,dzstar::Number) where {T}
  zr, zi = reim(z)
  drdx, didx = _dx_derivs(complex(dz),complex(dzstar))
  drdy, didy = _dy_derivs(complex(dz),complex(dzstar))
  return Dual{T}(zr,ForwardDiff.Partials((drdx,drdy))) +
      im*Dual{T}(zi,ForwardDiff.Partials((didx,didy)))
end

@inline ComplexComplexDual{T}(z::AbstractArray{S},dz::AbstractArray{S},dzstar::AbstractArray{S}) where {T,S<:Number} =
      map((u, v, w) -> ComplexComplexDual{T}(u,v,w),z,dz,dzstar)

@inline function ComplexRealDual{T}(z::Number,dz::Number) where {T}
  zr, zi = reim(z)
  dr, di = reim(dz)
  return Dual{T}(zr,dr) + im*Dual{T}(zi,di)
end

@inline ComplexRealDual{T}(z::AbstractArray{S},dz::AbstractArray{S}) where {T,S<:Number} =
      map((u, v) -> ComplexRealDual{T}(u,v),z,dz)

@inline ComplexComplexDual{T}() where {T} = ComplexComplexDual{T}(0.0,0.0,0.0)
@inline ComplexComplexDual{T}(z) where {T} = ComplexComplexDual{T}(z,0.0,0.0)
@inline ComplexComplexDual(args...) = ComplexComplexDual{Nothing}(args...)

@inline ComplexRealDual{T}() where {T} = ComplexRealDual{T}(0.0,0.0)
@inline ComplexRealDual{T}(z) where {T} = ComplexRealDual{T}(z,0.0)
@inline ComplexRealDual(args...) = ComplexRealDual{Nothing}(args...)

@inline Base.one(::Type{ComplexComplexDual{T}},z::Number) where T = ComplexComplexDual{T}(z,one(z),zero(z))
@inline Base.one(::Type{ComplexRealDual{T}},z::Number) where T = ComplexRealDual{T}(z,one(z))

@inline function value(z::ComplexDual)
    zr, zi = reim(z)
    return value(zr)+im*value(zi)
end

# extract the partial derivatives from a complex dual number
@inline function partials(z::ComplexDual)
    zr, zi = reim(z)
    return partials(zr), partials(zi)
end


# assemble the derivatives d/dz and d/dz* from the partials
@inline function _derivs(dvr::ForwardDiff.Partials{2},dvi::ForwardDiff.Partials{2})
    dvrx, dvry = dvr
    dvix, dviy = dvi
    return 0.5*(dvrx+dviy) + 0.5*im*(dvix-dvry), 0.5*(dvrx-dviy) + 0.5*im*(dvix+dvry)
end

# assemble the derivative from the partials
@inline function _derivs(dvr::ForwardDiff.Partials{1},dvi::ForwardDiff.Partials{1})
    return dvr[1] + im*dvi[1]
end

# get the d/dx derivatives from d/dz and d/dz*
@inline _dx_derivs(dvz::Complex{T},dvzstar::Complex{T}) where {T} = reim(   dvz +    dvzstar)

# get the d/dy derivatives from d/dz and d/dz*
@inline _dy_derivs(dvz::Complex{T},dvzstar::Complex{T}) where {T} = reim(im*dvz - im*dvzstar)




# preempt other function diffrules. The two entries correspond to d/dz and d/dz*
DiffRules.@define_diffrule Base.abs2(z) = :(conj($z)), :($z)
#DiffRules.@define_diffrule Base.conj(z) = :(0), :(1)
#DiffRules.@define_diffrule Base.real(z) = :(0.5), :(0.5)
#DiffRules.@define_diffrule Base.imag(z) = :(0.5im), :(-0.5im)

DiffRules.@define_diffrule Base.sqrt(z) = :(inv(2 * sqrt($z))), :(0)
DiffRules.@define_diffrule Base.log(z) = :(inv($z)), :(0)
DiffRules.@define_diffrule Base.:^(z,p) = :($p * ($z^($p - 1))), :(0)

DiffRules.@define_diffrule Base.abs(z) = :(0.5*conj($z)*inv(abs($z))), :(0.5*$z*inv(abs($z)))

# Extend functions of a single complex argument to accept dual argument
macro extend_unary_dual_to_complex(fcn)
    M = :Base
    f = :($fcn)
    dfz, dfzstar = DiffRules.diffrule(M,f,:v)
    fdef = quote
        function $M.$f(z::ComplexComplexDual{T}) where {T}
            v = value(z)
            dvz, dvzstar = extract_derivative(T,z)

            # chain rule
            dfdz =     $dfz*dvz     + $dfzstar*conj(dvzstar)
            dfdzstar = $dfz*dvzstar + $dfzstar*conj(dvz)

            return ComplexComplexDual{T}($M.$f(v),dfdz,dfdzstar)

        end
        function $M.$f(z::ComplexRealDual{T}) where {T}
            v = value(z)
            dv = extract_derivative(T,z)

            # chain rule
            df =     $dfz*dv     + $dfzstar*conj(dv)

            return ComplexRealDual{T}($M.$f(v),df)

        end
    end
    return esc(fdef)
end


# Extend functions of a single complex argument and another (non-dual) argument to accept
# dual argument in place of first argument
macro extend_binary_dual_to_complex(fcn)
    M = :Base
    f = :($fcn)
    dfz, dfzstar = DiffRules.diffrule(M,f,:v,:p)
    defs = quote end
    for R in AMBIGUOUS_TYPES
        expr = quote
            function $M.$f(z::ComplexComplexDual{T},p::$R) where {T}
                v = value(z)
                dvz, dvzstar = extract_derivative(T,z)

                # chain rule
                dfdz =     $dfz*dvz     + $dfzstar*conj(dvzstar)
                dfdzstar = $dfz*dvzstar + $dfzstar*conj(dvz)

                return ComplexComplexDual{T}($M.$f(v,p),dfdz,dfdzstar)

            end
            function $M.$f(z::ComplexRealDual{T},p::$R) where {T}
                v = value(z)
                dv = extract_derivative(T,z)

                # chain rule
                df =     $dfz*dv     + $dfzstar*conj(dv)

                return ComplexRealDual{T}($M.$f(v,p),df)

            end
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end


#@extend_unary_dual_to_complex conj
#@extend_unary_dual_to_complex real
#@extend_unary_dual_to_complex imag

@extend_unary_dual_to_complex abs2
@extend_unary_dual_to_complex sqrt
@extend_unary_dual_to_complex log
@extend_unary_dual_to_complex abs

@extend_binary_dual_to_complex ^

"""
    dualize(v::Vector{S},i::Int,T)

Return a Dual or complex Dual form of vector `v` (depending on whether `S` is `Float64`
or `ComplexF64`), with the partials of the `i`th component of the vector set to unit
values (i.e., to `1` or to `1,0`, respectively).
"""
function dualize(v::Vector{ComplexF64},i::Int,::Type{T}) where {T}
    @assert 1 <= i <= length(v) "Invalid index"
    d = one(ComplexComplexDual{T},v[i])
    dualv = convert(Vector{typeof(d)},v)
    dualv[i] = d
    return dualv
end

function dualize(v::Vector{Float64},i::Int,::Type{T}) where {T}
    @assert 1 <= i <= length(v) "Invalid index"
    d = Dual{T}(v[i],one(Float64))
    dualv = convert(Vector{typeof(d)},v)
    dualv[i] = d
    return dualv
end

"""
    ForwardDiff.derivative(f,z::Complex)

Compute the derivative of function `f` with respect to `z` and `conj(z)`
"""
@inline function derivative(f::F, z::C) where {F,C<:Complex}
    T = typeof(ForwardDiff.Tag(f, C))
    # making complex ensures that it gets dispatched to our extract_derivative
    # for complex Duals rather than the native one in ForwardDiff
    return extract_derivative(T,complex.(f(one(ComplexComplexDual{T},z))))
end

"""
    ForwardDiff.extract_derivative(T,d::ComplexDual)

Given a complex dual value `d` and tag `T`, extract the derivatives d/dz
and d/dz* from the partials in `d`.
"""
@inline function extract_derivative(::Type{T},d::ComplexDual{T}) where T
   dr, di = reim(d)
   return _derivs(partials(dr),partials(di))
end

# If the function outputs a tuple of duals, as would be the case for
# nested derivatives, then we deal with each one by one and assemble an overall
# list of derivatives
@inline function extract_derivative(::Type{T},dlist::NTuple{N,ComplexDual{T}}) where {N,T}
   out = ()
   for d in dlist
        out = (out...,extract_derivative(T,d)...)
    end
    out
end

# These are meant for array-valued outputs of functions, where each array
# element will be a derivative (or tuple of derivatives d/dz, d/dz*)
@inline extract_derivative(::Type{T},v::AbstractArray{<:ComplexComplexDual{T}}) where {T} =
        (d = map(x -> extract_derivative(T,x), v); return first.(d), last.(d))

@inline extract_derivative(::Type{T},v::AbstractArray{<:ComplexRealDual{T}}) where {T} =
        (d = map(x -> extract_derivative(T,x), v); return d)
