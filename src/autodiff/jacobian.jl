# ===== jacobian ===== #
import ForwardDiff: JACOBIAN_ERROR, reshape_jacobian


function jacobian(f, z::AbstractArray{Complex{V}},
            cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, z)) where {T, V}
    #CHK && checktag(T, f, z)
    checktag(T, f, z)
    vector_mode_jacobian(f, z, cfg)
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
            evalfcn=ForwardDiff.vector_mode_dual_eval) where {F,T,V,N}
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
    #dz_partials_wrap(ydual, nrange) = dz_partials(T, ydual, nrange)
    #dztmp, dzstartmp .= dz_partials_wrap.(ydual_reshaped, transpose(irange))
    for icol in irange, row in 1:size(dz, 1)
        dz[row, icol+offset], dzstar[row, icol+offset] = dz_partials(T, ydual[row], icol)
    end
    return dz, dzstar
end


function jacobian_chunk_mode_expr(work_array_definition::Expr, compute_ydual::Expr,
                                  dz_definition::Expr, dzstar_definition::Expr, y_definition::Expr)
    return quote
        @assert length(z) >= N "chunk size cannot be greater than length(z) ($(N) > $(length(z)))"

        # precalculate loop bounds
        zlen = length(z)
        remainder = zlen % N
        lastchunksize = ifelse(remainder == 0, N, remainder)
        lastchunkindex = zlen - lastchunksize + 1
        middlechunks = 2:div(zlen - lastchunksize, N)

        # seed work arrays
        $(work_array_definition)
        rseeds = cfg.rseeds
        iseeds = cfg.iseeds

        # do first chunk manually to calculate output type
        seed!(zdual, z, 1, rseeds, iseeds)
        $(compute_ydual)
        ydual isa AbstractArray || throw(JACOBIAN_ERROR)

        $(dz_definition)
        $(dzstar_definition)

        dz_reshaped = reshape_jacobian(dz, ydual, zdual)
        dzstar_reshaped = reshape_jacobian(dzstar, ydual, zdual)

        extract_jacobian_chunk!(T, dz_reshaped, dzstar_reshaped, ydual, 1, N)
        seed!(zdual, z, 1)

        # do middle chunks
        for c in middlechunks
            i = ((c - 1) * N + 1)
            seed!(zdual, z, i, rseeds,iseeds)
            $(compute_ydual)
            extract_jacobian_chunk!(T, dz_reshaped, dzstar_reshaped, ydual, i, N)
            seed!(zdual, z, i)
        end

        # do final chunk
        seed!(zdual, z, lastchunkindex, rseeds, iseeds, lastchunksize)
        $(compute_ydual)
        extract_jacobian_chunk!(T, dz_reshaped, dzstar_reshaped, ydual, lastchunkindex, lastchunksize)

        $(y_definition)

        return dz, dzstar
    end
end

@eval function chunk_mode_jacobian(f::F, z, cfg::ComplexGradientConfig{T,V,N}) where {F,T,V,N}
    $(jacobian_chunk_mode_expr(quote
                                   zdual = cfg.duals
                                   seed!(zdual, z)
                               end,
                               :(ydual = f(zdual)),
                               :(dz = similar(z, Complex{valtype(ydual)},length(ydual), zlen)),
                               :(dzstar = similar(z, Complex{valtype(ydual)},length(ydual), zlen)),
                               :()))
end
