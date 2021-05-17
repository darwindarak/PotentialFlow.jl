
"""
    ForwardDiff.derivative(f,z::Complex)

Compute the derivative of function `f` with respect to `z` and `conj(z)`
"""
@inline function derivative(f::F, z::C) where {F,C<:Complex}
    T = typeof(ForwardDiff.Tag(f, C))
    # making f output complex ensures that it gets dispatched to our extract_derivative
    # for complex Duals rather than the native one in ForwardDiff
    return extract_derivative(T,complex.(f(one(ComplexDual{T},z))))
end

"""
    ForwardDiff.extract_derivative(T,d::ComplexDual)

Given a complex dual value `d` and tag `T`, extract the derivatives d/dz
and d/dz* from the partials in `d`.
"""
@inline function extract_derivative(::Type{T},d::ComplexDual{T}) where T
   return dz_partials(d,1)
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
