
using ForwardDiff
using DiffRules

import ForwardDiff:value,partials,derivative,extract_derivative
const AMBIGUOUS_TYPES = (AbstractFloat, Irrational, Integer, Rational, Real, RoundingMode, ComplexF64)

# preempt other function diffrules. The two entries correspond to d/dz and d/dz*
DiffRules.@define_diffrule Base.abs2(z) = :(conj($z)), :($z)
DiffRules.@define_diffrule Base.conj(z) = :(0), :(1)
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
        function $M.$f(z::Complex{<:ForwardDiff.Dual{T}}) where {T}
            zr, zi = reim(z)
            #v = value(zr)+im*value(zi)
            v = value(z)
            dvz, dvzstar = _dz_derivs(partials(zr),partials(zi))

            # chain rule
            dfdz =     $dfz*dvz     + $dfzstar*conj(dvzstar)
            dfdzstar = $dfz*dvzstar + $dfzstar*conj(dvz)

            return complex_dual(T,$M.$f(v),dfdz,dfdzstar)
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
            function $M.$f(z::Complex{<:ForwardDiff.Dual{T}},p::$R) where {T}
                zr, zi = reim(z)
                #v = value(zr)+im*value(zi)
                v = value(z)
                dvz, dvzstar = _dz_derivs(partials(zr),partials(zi))

                # chain rule
                dfdz =     $dfz*dvz     + $dfzstar*conj(dvzstar)
                dfdzstar = $dfz*dvzstar + $dfzstar*conj(dvz)

                return complex_dual(T,$M.$f(v,p),dfdz,dfdzstar)
            end
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end

@inline function _dz_derivs(dvr::ForwardDiff.Partials{2},dvi::ForwardDiff.Partials{2})
    dvrx, dvry = dvr
    dvix, dviy = dvi
    return 0.5*(dvrx+dviy) + 0.5*im*(dvix-dvry), 0.5*(dvrx-dviy) + 0.5*im*(dvix+dvry)
end

@inline _dx_derivs(dvz::Complex{T},dvzstar::Complex{T}) where {T} = reim(   dvz +    dvzstar)

@inline _dy_derivs(dvz::Complex{T},dvzstar::Complex{T}) where {T} = reim(im*dvz - im*dvzstar)

@inline function value(z::Complex{<:ForwardDiff.Dual{T}}) where {T}
    zr, zi = reim(z)
    return value(zr)+im*value(zi)
end

# extract the partial derivatives from a complex dual number
@inline function partials(z::Complex{<:ForwardDiff.Dual{T}}) where {T}
    zr, zi = reim(z)
    return partials(zr), partials(zi)
end

@inline function complex_dual(::Type{T},z,dz,dzstar) where {T}
    zr, zi = reim(z)
    drdx, didx = _dx_derivs(dz,dzstar)
    drdy, didy = _dy_derivs(dz,dzstar)
    return ForwardDiff.Dual{T}(zr,ForwardDiff.Partials((drdx,drdy))) +
        im*ForwardDiff.Dual{T}(zi,ForwardDiff.Partials((didx,didy)))
end

@inline complex_dual(d::Complex{<:ForwardDiff.Dual{T}}) where {T} = d


@extend_unary_dual_to_complex conj
#@extend_unary_dual_to_complex real
#@extend_unary_dual_to_complex imag
@extend_unary_dual_to_complex abs2
@extend_unary_dual_to_complex sqrt
@extend_unary_dual_to_complex log
@extend_unary_dual_to_complex abs

@extend_binary_dual_to_complex ^


"""
    ForwardDiff.derivative(f,z::Complex)

Compute the derivative of function `f` with respect to `z` and `conj(z)`
"""
@inline function ForwardDiff.derivative(f::F, z::C) where {F,C<:Complex}
    T = typeof(ForwardDiff.Tag(f, C))
    outr, outi = reim(f(complex_dual(T,z,one(z),zero(z))))
    return _dz_derivs(partials(outr),partials(outi))
end
