# ===== jacobian ===== #
import ForwardDiff: JACOBIAN_ERROR, reshape_jacobian, chunksize


function jacobian(f, z::AbstractArray{Complex{V}},
            cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, z)) where {T, V}
    #CHK && checktag(T, f, z)
    checktag(T, f, z)
    if chunksize(cfg) == length(z)
      vector_mode_jacobian(f, z, cfg)
    else
      chunk_mode_jacobian(f, z, cfg)
    end
end

function extract_jacobian!(::Type{T}, dz::AbstractArray, dzstar::AbstractArray,
                ydual::AbstractArray{<:ComplexDual{T,V,N}}, n) where {T,V,N}
    @assert size(dz) == size(dzstar) "Inconsistent result sizes"

    dz_reshaped = reshape(dz, length(ydual), n)
    dzstar_reshaped = reshape(dzstar, length(ydual), n)
    for col in 1:size(dz_reshaped, 2), row in 1:size(dz_reshaped, 1)
        dz_reshaped[row, col], dzstar_reshaped[row, col] = dz_partials(T, ydual[row], col)
    end
    return dz, dzstar
end



function vector_mode_jacobian(f::F, z, cfg::ComplexGradientConfig{T,V,N},
            evalfcn=ForwardDiff.vector_mode_dual_eval!) where {F,T,V,N}
    ydual = evalfcn(f, z, cfg)
    dz = similar(ydual, Complex{valtype(eltype(ydual))}, length(ydual), N)
    dzstar = similar(ydual, Complex{valtype(eltype(ydual))}, length(ydual), N)
    extract_jacobian!(T, dz, dzstar, ydual, N)
    return dz, dzstar
end


function extract_jacobian_chunk!(::Type{T}, dz, dzstar, ydual, index, chunksize) where {T}
    ydual_reshaped = vec(ydual)
    offset = index - 1
    irange = 1:chunksize
    col = irange .+ offset

    # Use closure to avoid GPU broadcasting with Type
    dz_partials_wrap(ydual, nrange) = dz_partials(T, ydual, nrange)
    dztmp = collect(Iterators.flatten(dz_partials_wrap.(ydual_reshaped, transpose(irange))))
    len, ylen = length(dztmp), length(ydual_reshaped)

    dz[:,col] .= reshape(view(dztmp,1:2:len-1),ylen,chunksize)
    dzstar[:,col] .= reshape(view(dztmp,2:2:len),ylen,chunksize)

    #for icol in irange, row in 1:size(dz, 1)
    #    dz[row, icol+offset], dzstar[row, icol+offset] = dz_partials(T, ydual[row], icol)
    #end
    return dz, dzstar
end


@eval function chunk_mode_jacobian(f::F, z, cfg::ComplexGradientConfig{T,V,N}) where {F,T,V,N}
    $(chunk_mode_gradient_expr(quote
                                   zdual = cfg.duals
                                   seed!(zdual, z)
                               end,
                               quote
                                 dz = similar(z, Complex{valtype(ydual)},length(ydual), zlen)
                                 dzstar = similar(z, Complex{valtype(ydual)},length(ydual), zlen)
                                 dz_reshaped = reshape_jacobian(dz, ydual, zdual)
                                 dzstar_reshaped = reshape_jacobian(dzstar, ydual, zdual)
                               end,
                               :(dz, dzstar),
                               :(ydual = f(zdual)),
                               :(ydual isa AbstractArray || throw(JACOBIAN_ERROR)),
                               :(seed!(zdual, z, index, rseeds, iseeds, chunksize)),
                               :(seed!(zdual, z, index)),
                               :(extract_jacobian_chunk!(T, dz_reshaped, dzstar_reshaped, ydual, index, chunksize)),
                               :()))
end
