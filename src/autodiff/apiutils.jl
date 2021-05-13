@generated function construct_complex_seeds(::Type{Partials{N,V}}) where {N,V}
    return Expr(:tuple, Expr(:tuple,[:(single_seed(Partials{N,V}, Val{2*$i-1}())) for i in 1:N÷2]...),
                        Expr(:tuple,[:(single_seed(Partials{N,V}, Val{2*$i-0}())) for i in 1:N÷2]...))
end

function seed!(duals::AbstractArray{<:ComplexDual{T,V,M}}, x,
               rseed::Partials{M,V} = zero(Partials{M,V}),iseed::Partials{M,V} = zero(Partials{M,V})) where {T,V,N,M}
    for i in eachindex(duals)
        duals[i] = ComplexDual{T,V,M}(x[i],rseed,iseed)
    end
    return duals
end

function seed!(duals::AbstractArray{<:ComplexDual{T,V,M}}, x,
               rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}}) where {T,V,N,M}
    for i in 1:N
        duals[i] = ComplexDual{T,V,M}(x[i],rseeds[i],iseeds[i])
    end
    return duals
end

function seed!(duals::AbstractArray{<:ComplexDual{T,V,M}}, x, index,
               rseed::Partials{M,V} = zero(Partials{M,V}),iseed::Partials{M,V} = zero(Partials{M,V})) where {T,V,N,M}
    offset = index - 1
    dual_inds = (1:M÷2) .+ offset
    duals[dual_inds] .= ComplexDual{T,V,M}.(view(x, dual_inds), Ref(rseed),Ref(iseed))
    return duals
end

function seed!(duals::AbstractArray{<:ComplexDual{T,V,M}}, x, index,
               rseeds::NTuple{N,Partials{M,V}},iseeds::NTuple{N,Partials{M,V}}, chunksize = N) where {T,V,N,M}
    offset = index - 1
    seed_inds = 1:chunksize
    dual_inds = seed_inds .+ offset
    duals[dual_inds] .= ComplexDual{T,V,M}.(view(x, dual_inds), getindex.(Ref(rseeds), seed_inds), getindex.(Ref(iseeds), seed_inds))
    return duals
end



function dualize(::Type{T},x::AbstractArray{<:Complex{V}}) where {T,V}
    xdual = similar(x,ComplexDual{T,V,2*length(v)})
    seed!(xdual,x)
end

# =====  other seed functions ==== ##

"""
    seed(v::Number,cfg::ComplexGradientConfig)

Given a number `v`, create a dual-valued copy of `v` with
a unit seed. The use of complex duals (even if `v` is real-valued)
ensures correct treatment for gradient calculation, since some
calculations will be complex-valued.
"""
function seed(v::Complex,cfg::ComplexGradientConfig)
    vduals = copy(cfg.duals)
    seed!(vduals,[v],cfg.rseeds,cfg.iseeds)
    return vduals[1]
end

function seed(v::Real,cfg::ComplexGradientConfig)
    vduals = copy(cfg.duals)
    seed!(vduals,complex([v]),cfg.rseeds,cfg.iseeds)
    return real(vduals[1])
end
