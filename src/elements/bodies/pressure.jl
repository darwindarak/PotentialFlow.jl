"""
    pressure(ζ,v::Vector{Element},b::ConformalBody)

Return the pressure at ζ (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`.
"""
function pressure(ζ,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
        out = real(zero(ζ))

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇
        Ũ̇, Ṽ̇ = reim(b.c̈*exp(-im*b.α))
        Ω̇ = b.α̈

        for (j,vj) in enumerate(v)
            zj,Γj  = Elements.position(vj), circulation(vj)
            out -= 0.5*Γj^2*Πvv(ζ,vj,vj,b;kwargs...)
            for vk in v[1:j-1]
                zk,Γk  = Elements.position(vk), circulation(vk)
                out -= Γj*Γk*Πvv(ζ,vj,vk,b;kwargs...)
            end

            out -= Γj*Ũ*ΠUv(ζ,vj,b;kwargs...)
            out -= Γj*Ṽ*ΠVv(ζ,vj,b;kwargs...)
            out -= Γj*Ω*ΠΩv(ζ,vj,b;kwargs...)
        end

        out -= 0.5*Ũ^2*ΠUU(ζ,b)
        out -= 0.5*Ṽ^2*ΠVV(ζ,b)
        out -= 0.5*Ω^2*ΠΩΩ(ζ,b)
        out -= Ũ*Ṽ*ΠUV(ζ,b)
        out -= Ũ*Ω*ΠUΩ(ζ,b)
        out -= Ṽ*Ω*ΠVΩ(ζ,b)

        out -= Ũ̇*PU(ζ,b)
        out -= Ṽ̇*PV(ζ,b)
        out -= Ω̇*PΩ(ζ,b)


        return out
end

"""
    dphivdzv(ζ,b::ConformalBody,v::Element)

Calculates the derivative of the scalar potential field ``\\phi`` with respect to a change in
the body-fixed coordinates of element `v`, also accounting for the change in the image of `v`.
"""
function dphivdzv(ζ,b::Bodies.ConformalBody,v::Element)
    return conj(Bodies.transform_velocity(_dphivdzetav_conj(ζ,b,v),v.z,b))
end

function _dphivdzetav_conj(ζ,b::Bodies.ConformalBody,v::Element)
    out = induce_velocity(ζ,_unit_strength_copy(v),0.0)
    b_img = _set_image_for_stationary_body(b,v,Val(false))
    out -= conj(induce_velocity(ζ,b_img,0.0)/Elements.position(v)^2)
    return -0.5*out
end

"""
    wvv(targ::Element,src::Element,b::ConformalBody)

Calculate the complex velocity ``w = u - iv`` of the target element `targ` due to element `src` with unit strength,
also accounting for the images of `src`. If `targ` and `src` are the same element, then
it makes the Routh correction.
"""
@inline wvv(targ::Element,src::Element,b::Bodies.ConformalBody) = wv(targ.z,b,src) +
                                                                   conj(Bodies.transform_velocity(_routh_conj(targ.z,b,Val(targ==src)),targ.z,b))

function _routh_conj(ζ,b,::Val{true})
    dz̃, ddz̃ = b.dm(ζ)
    return -conj(ddz̃/(4π*im*dz̃))
end

_routh_conj(ζ,b,::Val{false}) = complex(0.0)

#=
function P(ζ,v::Element,b::Bodies.ConformalBody)
    return abs2.(wv(ζ,b,v)) + 4.0*real(dphivdzv(ζ,b,v)*conj(wvv(v,v,b)))
end
=#

function Πvv(ζ,targ::Element,src::Element,b::Bodies.ConformalBody)
    return real(wv(ζ,b,targ).*conj(wv(ζ,b,src))) +
           2.0*real(dphivdzv(ζ,b,targ)*conj(wvv(targ,src,b))) +
           2.0*real(dphivdzv(ζ,b,src)*conj(wvv(src,targ,b)))
end

function ΠUv(ζ,v::Element,b::Bodies.ConformalBody)
    return -real(winf1(ζ,b).*conj(wv(ζ,b,v))) -
           4.0*real(dphivdzv(ζ,b,v)*conj(winf1(v.z,b)))
end

function ΠVv(ζ,v::Element,b::Bodies.ConformalBody)
    return -real(winf2(ζ,b).*conj(wv(ζ,b,v))) -
           4.0*real(dphivdzv(ζ,b,v)*conj(winf2(v.z,b)))
end

function ΠΩv(ζ,v::Element,b::Bodies.ConformalBody)
    return -real(wrinf(ζ,b).*conj(wv(ζ,b,v))) -
           4.0*real(dphivdzv(ζ,b,v)*conj(wrinf(v.z,b)))
end

ΠUU(ζ,b::Bodies.ConformalBody) = abs2.(winf1(ζ,b)) .- 1.0
ΠVV(ζ,b::Bodies.ConformalBody) = abs2.(winf2(ζ,b)) .- 1.0
ΠΩΩ(ζ,b::Bodies.ConformalBody) = abs2.(wrinf(ζ,b)) .- abs2.(b.m.(ζ))

ΠUV(ζ,b::Bodies.ConformalBody) = real(winf1(ζ,b).*conj(winf2(ζ,b)))
ΠUΩ(ζ,b::Bodies.ConformalBody) = real(winf1(ζ,b).*conj(wrinf(ζ,b))) + imag(b.m.(ζ))
ΠVΩ(ζ,b::Bodies.ConformalBody) = real(winf2(ζ,b).*conj(wrinf(ζ,b))) - real(b.m.(ζ))

PU(ζ,b::Bodies.ConformalBody) = real(Fbt1(ζ,b))
PV(ζ,b::Bodies.ConformalBody) = real(Fbt2(ζ,b))
PΩ(ζ,b::Bodies.ConformalBody) = real(Fbr(ζ,b))
