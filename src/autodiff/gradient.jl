
# Given a function f of z and a complex array z, evaluate the gradient wrt z at z
function gradient(f, z::AbstractArray{Complex{V}},
            cfg::ComplexGradientConfig{T} = ComplexGradientConfig(f, z)) where {T, V}
    #CHK && checktag(T, f, z)
    checktag(T, f, z)
    vector_mode_gradient(f, z, cfg)
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

function extract_gradient_chunk!(::Type{T}, dz, dzstar, d::ComplexDual{T,V,N}, index, chunksize) where {T,V,N}
    offset = index - 1
    for i in 1:chunksize
        tmp, tmpstar = dz_partials(T, d, i)
        dz[i + offset] = tmp
        dzstar[i + offset] = tmpstar
    end
    return dz, dzstar
end


function vector_mode_gradient(f::F, z, cfg::ComplexGradientConfig{T},
                evalfcn=ForwardDiff.vector_mode_dual_eval) where {F,T}
    ydual = evalfcn(f, z, cfg)
    dz = similar(z, Complex{valtype(ydual)})
    dzstar = similar(z, Complex{valtype(ydual)})
    return extract_gradient!(T, dz, dzstar,ydual)
end

function ForwardDiff.vector_mode_dual_eval(f::F, z, cfg::ComplexGradientConfig) where {F}
    zduals = cfg.duals
    seed!(zduals,z,cfg.rseeds,cfg.iseeds)
    return f(zduals)
end

function chunk_mode_complex_gradient_expr(dz_definition::Expr, dzstar_definition::Expr)
  return quote

    @assert length(z) >= N "chunk size cannot be greater than length(z) ($(N) > $(length(z)))"

    # precalculate loop bounds
    zlen = length(z)
    remainder = zlen % N
    lastchunksize = ifelse(remainder == 0, N, remainder)
    lastchunkindex = zlen - lastchunksize + 1
    middlechunks = 2:div(zlen - lastchunksize, N)

    # seed work vectors
    zdual = cfg.duals
    rseeds = cfg.rseeds
    iseeds = cfg.iseeds
    seed!(zdual, z)
    seed!(zdual, z, 1, rseeds,iseeds)
    ydual = f(zdual)
    $(dz_definition)
    $(dzstar_definition)

    extract_gradient_chunk!(T,dz,dzstar,ydual,1,N)
    seed!(zdual, z, 1)

    for c in middlechunks
      i = ((c - 1) * N + 1)
      seed!(zdual, z, i, rseeds,iseeds)
      ydual = f(zdual)
      extract_gradient_chunk!(T, dz, dzstar, ydual, i, N)
      seed!(zdual, z, i)
    end

    seed!(zdual, z, lastchunkindex, rseeds, iseeds, lastchunksize)
    ydual = f(zdual)
    extract_gradient_chunk!(T, dz, dzstar, ydual, lastchunkindex, lastchunksize)

    return dz, dzstar
  end
end

@eval function chunk_mode_gradient(f::F, z, cfg::ComplexGradientConfig{T,V,N}) where {F,T,V,N}
    $(chunk_mode_complex_gradient_expr(:(dz = similar(z, Complex{valtype(ydual)})),
                                       :(dzstar = similar(z, Complex{valtype(ydual)}))))
end
