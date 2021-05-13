# =====  Configurations =====  #
"""
    ComplexGradientConfig(f,z)

Create a configuration for calculation of the gradient of a function
`f` with respect to complex array `z` by automatic differentiation.
"""
struct ComplexGradientConfig{T,V,N,D,M} <: ForwardDiff.AbstractConfig{N}
    rseeds::NTuple{N,Partials{M,V}}
    iseeds::NTuple{N,Partials{M,V}}
    duals::D
end

function ComplexGradientConfig(f::F,
                        x::AbstractArray{Complex{V}},
                        ::Chunk{N} = Chunk(x,length(x)),
                        ::T = Tag(f, Complex{V})) where {F,V,N,T}
    rseeds, iseeds = construct_complex_seeds(Partials{2N,V})
    duals = similar(x, ComplexDual{T,V,2N})
    return ComplexGradientConfig{T,V,N,typeof(duals),2N}(rseeds, iseeds, duals)
end
