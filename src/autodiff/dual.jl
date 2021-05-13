# ComplexDual encompasses all complex dual types
const ComplexDual{T,V,N} = Complex{Dual{T,V,N}}

# Constructor for ComplexDual based on real Partials for dreal and dimag
@inline function ComplexDual{T,V,N}(z::Number,dr::Partials{N,V},di::Partials{N,V}) where {T,N,V<:Real}
    zr, zi = reim(z)
    return Dual{T}(zr,dr) + im*Dual{T}(zi,di)
end

@inline ComplexDual{T}(z::Number,u::Partials{N,V},v::Partials{N,V}) where {T,N,V} =
    ComplexDual{T,V,N}(z,u,v)

# Constructor for ComplexDual based on complex Partials for dz and dz*
@inline function ComplexDual{T}(z::Number,dz::Partials{N,V},dzstar::Partials{N,V}) where {T,N,V<:Complex}
    dr, di = reim_partials(dz,dzstar)
    return ComplexDual{T}(z,dr,di)
end

@inline ComplexDual{T}(z::Number,dz::Complex,dzstar::Complex) where {T} =
    ComplexDual{T}(z,Partials(tuple(dz)),Partials(tuple(dzstar)))

@inline Base.one(::Type{<:ComplexDual{T}},z::Number) where T = ComplexDual{T}(z,one(z),zero(z))


@inline valtype(::ComplexDual{T,V,N}) where {T,V,N} = V
@inline valtype(::Type{ComplexDual{T,V,N}}) where {T,V,N} = V
@inline valtype(::Array{<:ComplexDual{T,V,N}}) where {T,V,N} = V
@inline valtype(::Type{<:Array{<:ComplexDual{T,V,N}}}) where {T,V,N} = V


# ==== extensions of Partials ==== #
@inline Base.conj(v::Partials{N,V}) where {N,V} = Partials{N,V}(conj.(v.values))
@inline Base.real(v::Partials{N,Complex{V}}) where {N,V} = Partials{N,V}(real.(v.values))
@inline Base.imag(v::Partials{N,Complex{V}}) where {N,V} = Partials{N,V}(imag.(v.values))

@inline Base.:*(x::Complex, partials::Partials) = partials*x

@inline function Base.:*(partials::Partials, x::Complex)
        return Partials(ForwardDiff.scale_tuple(partials.values, x))
end

# extract the partial derivatives from a complex dual number
@inline function ForwardDiff.partials(d::ComplexDual)
    dr, di = reim(d)
    return partials(dr), partials(di)
end

# assemble the derivatives d/dz and d/dz* from the partials
@inline function _derivs(dvrx,dvry,dvix,dviy)
    return 0.5*(dvrx+dviy) + 0.5*im*(dvix-dvry), 0.5*(dvrx-dviy) + 0.5*im*(dvix+dvry)
end

# Given a complex dual number, extract Partials of dz and dz*
@inline function dz_partials(d::ComplexDual{T,V,N}) where {T,V,N}
    dr, di = partials(d)
    tmp  = [_derivs(dr[i], dr[i+1], di[i],di[i+1]) for i in 1:2:N]
    return Partials(tuple(first.(tmp)...)), Partials(tuple(last.(tmp)...))
end

# For the dz and dz* partials with respect to the ith element
@inline Base.@propagate_inbounds function dz_partials(d::ComplexDual, i)
    dr, di = partials(d)
    return _derivs(dr[2i-1], dr[2i], di[2i-1],di[2i])
end

@inline Base.@propagate_inbounds dz_partials(::Type{T}, d::ComplexDual{T}, i...) where T = dz_partials(d, i...)

# For dealing with arrays of complex duals, return the dz and dz*
# partials as separate arrays
@inline dz_partials(d::AbstractVector{ComplexDual{T,V,N}}) where {T,V,N} =
  (tmp = dz_partials.(d); return first.(tmp), last.(tmp))

# For the dz and dz* partials with respect to the ith element
@inline Base.@propagate_inbounds dz_partials(d::AbstractVector{ComplexDual{T,V,N}},i) where {T,V,N} =
  (tmp = dz_partials.(d,i); return first.(tmp), last.(tmp))

@inline Base.@propagate_inbounds dz_partials(::Type{T}, d::AbstractVector{ComplexDual{T,V,N}}, i...) where {T,V,N} = dz_partials(d, i...)

# Given complex Partials for d/dz and d/dz*, return partials for dr/dx,y and di/dx,y
@inline function reim_partials(dz::Partials{N,V},dzstar::Partials{N,V}) where {T,V<:Complex,N}
    drx, dix = reim(dz+dzstar)
    dry, diy = reim(im*dz-im*dzstar)
    return _splice_xy(drx,dry), _splice_xy(dix,diy)
end

# Given partials for the Create separate partials for the real part and the imaginary part.
# These are organized as drx[1],dry[1],drx[2],dry[2],... and dix[1],diy[1],dix[2],diy[2],...
@inline function _splice_xy(dx::Partials{N},dy::Partials{N}) where {N}
    return Partials(tuple([ifelse(mod(i,2) == 0,dy[max(i÷2,1)],dx[(i+1)÷2]) for i in 1:2N]...))
end

# Return the value of a complex dual
@inline function value(d::ComplexDual)
    dr, di = reim(d)
    return value(dr)+im*value(di)
end
