"""
    pressure(ζ,v::Vector{Element},b::ConformalBody)

Return the pressure at ζ (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`.
"""
function pressure(ζ,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
        out = real(zero(ζ))
        for (j,vj) in enumerate(v)
            zj,Γj  = Elements.position(vj), circulation(vj)
            out -= 0.5*Γj^2*P(ζ,vj,b;kwargs...)
            for vk in v[1:j-1]
                zk,Γk  = Elements.position(vk), circulation(vk)
                out -= Γj*Γk*Π(ζ,vj,vk,b;kwargs...)
            end
        end
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

function P(ζ,v::Element,b::Bodies.ConformalBody)
    return abs2.(wv(ζ,b,v)) + 4.0*real(dphivdzv(ζ,b,v)*conj(wvv(v,v,b)))
end

function Π(ζ,targ::Element,src::Element,b::Bodies.ConformalBody)
    return real(wv(ζ,b,targ).*conj(wv(ζ,b,src))) +
           2.0*real(dphivdzv(ζ,b,targ)*conj(wvv(targ,src,b))) +
           2.0*real(dphivdzv(ζ,b,src)*conj(wvv(src,targ,b)))
end
