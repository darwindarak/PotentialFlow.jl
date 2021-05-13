# ===== jacobian ===== #

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
