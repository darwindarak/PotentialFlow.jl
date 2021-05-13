# ======== Extension of diff rules to complex functions ==========#


# preempt other function diffrules. The two entries correspond to d/dz and d/dz*
DiffRules.@define_diffrule Base.abs2(z) = :(conj($z)), :($z)
DiffRules.@define_diffrule Base.sqrt(z) = :(inv(2 * sqrt($z))), :(0)
DiffRules.@define_diffrule Base.log(z) = :(inv($z)), :(0)
DiffRules.@define_diffrule Base.:^(z,p) = :($p * ($z^($p - 1))), :(0)

DiffRules.@define_diffrule Base.abs(z) = :(0.5*conj($z)*inv(abs($z))), :(0.5*$z*inv(abs($z)))

# ======== Train functions to deal with complex duals ==========#


# Extend functions of a single complex argument to accept dual argument
macro extend_unary_dual_to_complex(fcn)
    M = :Base
    f = :($fcn)
    dfz, dfzstar = DiffRules.diffrule(M,f,:v)
    fdef = quote
        function $M.$f(z::ComplexDual{T}) where {T}
            v = value(z)
            dvz, dvzstar = dz_partials(z)

            # chain rule
            dfdz =     $dfz*dvz     + $dfzstar*conj(dvzstar)
            dfdzstar = $dfz*dvzstar + $dfzstar*conj(dvz)

            return ComplexDual{T}($M.$f(v),dfdz,dfdzstar)
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
            function $M.$f(z::ComplexDual{T},p::$R) where {T}
                v = value(z)
                dvz, dvzstar = dz_partials(z)

                # chain rule
                dfdz =     $dfz*dvz     + $dfzstar*conj(dvzstar)
                dfdzstar = $dfz*dvzstar + $dfzstar*conj(dvz)

                return ComplexDual{T}($M.$f(v,p),dfdz,dfdzstar)

            end
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end

@extend_unary_dual_to_complex abs2
@extend_unary_dual_to_complex sqrt
@extend_unary_dual_to_complex log
@extend_unary_dual_to_complex abs

@extend_binary_dual_to_complex ^
