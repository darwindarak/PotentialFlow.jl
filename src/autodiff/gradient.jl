
# Given a function f of z and a complex array z, evaluate the gradient wrt z at z
function gradient(f, z::AbstractArray{Complex{V}},
            cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, z)) where {T, V}
    #CHK && checktag(T, f, z)
    checktag(T, f, z)
    if chunksize(cfg) == length(z)
      vector_mode_gradient(f, z, cfg)
    else
      chunk_mode_gradient(f, z, cfg)
    end
end

# Given a complex dual, return dz and dzstar as Partial types. Ugly way of doing
# this
@inline function extract_gradient!(dz,dzstar,d::ComplexDual{T,V,N}) where {T,V,N}
    tmp, tmpstar = dz_partials(d)
    dz .= tmp
    dzstar .= tmpstar
    return dz,dzstar
end

# This is for making sure that the desired gradient has the same tag as the dual
@inline extract_gradient!(::Type{T},dz,dzstar,d::ComplexDual{T,V,N}) where {T,V,N} =
          extract_gradient!(dz,dzstar,d)


function vector_mode_gradient(f::F, z, cfg::ComplexGradientConfig{T},
                evalfcn=ForwardDiff.vector_mode_dual_eval!) where {F,T}
    ydual = evalfcn(f, z, cfg)
    dz = similar(z, Complex{valtype(ydual)})
    dzstar = similar(z, Complex{valtype(ydual)})
    return extract_gradient!(T, dz, dzstar,ydual)
end

function extract_gradient_chunk!(::Type{T}, dz, dzstar, d::ComplexDual{T,V,N}, index, chunksize) where {T,V,N}
    offset = index - 1
    for i in 1:chunksize
        tmp, tmpstar = dz_partials(T, d, i)
        dz[i + offset] = tmp
        dzstar[i + offset] = tmpstar
    end
    return dz, dzstar
end




function chunk_mode_gradient_expr(work_array_definition::Expr,result_definition::Expr,
                                          result::Expr,compute_ydual::Expr,
                                          error_definition::Expr,
                                          seed_chunk_definition::Expr, unseed_definition::Expr,
                                          extract_definition::Expr,y_definition::Expr)
  return quote

    @assert length(z) >= N "chunk size cannot be greater than length(z) ($(N) > $(length(z)))"

    # precalculate loop bounds
    zlen = length(z)
    remainder = zlen % N
    lastchunksize = ifelse(remainder == 0, N, remainder)
    lastchunkindex = zlen - lastchunksize + 1
    middlechunks = 2:div(zlen - lastchunksize, N)

    # seed work vectors
    $(work_array_definition)
    rseeds = cfg.rseeds
    iseeds = cfg.iseeds

    index, chunksize = 1, N
    $(seed_chunk_definition)
    $(compute_ydual)
    $(error_definition)
    $(result_definition)

    $(extract_definition)
    $(unseed_definition)

    for c in middlechunks
      i = ((c - 1) * N + 1)
      index, chunksize = i, N
      $(seed_chunk_definition)
      $(compute_ydual)
      $(extract_definition)
      $(unseed_definition)
    end
    index, chunksize = lastchunkindex, lastchunksize
    $(seed_chunk_definition)
    $(compute_ydual)
    $(extract_definition)

    $(y_definition)

    return $(result)
  end
end

@eval function chunk_mode_gradient(f::F, z, cfg::ComplexGradientConfig{T,V,N}) where {F,T,V,N}
    $(chunk_mode_gradient_expr(quote
                                  zdual = cfg.duals
                                  seed!(zdual, z)
                                end,
                                quote
                                  dz = similar(z, Complex{valtype(ydual)})
                                  dzstar = similar(z, Complex{valtype(ydual)})
                                end,
                                :(dz, dzstar),
                                :(ydual = f(zdual)),
                                :(),
                                :(seed!(zdual, z, index, rseeds, iseeds, chunksize)),
                                :(seed!(zdual, z, index)),
                                :(extract_gradient_chunk!(T, dz, dzstar, ydual, index, chunksize)),
                                :()))
end
