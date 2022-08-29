#=
 Routines that provide the unit basis fields for a given body
 Note that the unit velocity fields are provided in standard
 complex form, w = u - iv.

 The conventions of these routines:
 wc to denote ŵ (velocity in the circle plane)
 ζ to renote a point in the circle plane
=#

export Fbt1, Fbt2, Fbr, Fv, wbt1, wbt2, wbr, winf1, winf2, wrinf, wv,
        wcbt1, wcbt2, wcbr, wcv

#### Basis complex potentials ####

"""
    Fbt1(ζ,b::ConformalBody)

Calculate the complex potential at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline Fbt1(ζ,b::ConformalBody) = Elements.complexpotential(ζ,_set_motion_in_quiescent(b,1.0,0.0,0.0))

"""
    Fbt2(ζ,b::ConformalBody)

Calculate the complex potential at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the ỹ direction (of its own coordinate system).
"""
@inline Fbt2(ζ,b::ConformalBody) = Elements.complexpotential(ζ,_set_motion_in_quiescent(b,0.0,1.0,0.0))

"""
    Fbr(ζ,b::ConformalBody)

Calculate the complex potential at `ζ` in the circle plane due to motion of body `b` with unit angular velocity
about its reference point
"""
@inline Fbr(ζ,b::ConformalBody) = Elements.complexpotential(ζ,_set_motion_in_quiescent(b,0.0,0.0,1.0))

"""
    Fv(ζ,b::ConformalBody,v::Element)

Calculate the complex potential at `ζ` in the circle plane due to a unit-strength vortex at the location of `v`,
along with its images inside the body `b`. It is presumed that the location of `v` is given in the circle plane.
"""
@inline Fv(ζ,b::ConformalBody,v::Element) = Elements.complexpotential(ζ,(_set_image_for_stationary_body(b,v,Val(true)),_unit_strength_copy(v)))

#### Basis velocities in circle plane ####

"""
    wcbt1(ζ,b::ConformalBody)

Calculate the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline wcbt1(ζ,b::ConformalBody) = conj(_wcbt1_conj(ζ,b))

"""
    wcbt2(ζ,b::ConformalBody)

Calculate the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline wcbt2(ζ,b::ConformalBody) = conj(_wcbt2_conj(ζ,b))

"""
    wcbr(ζ,b::ConformalBody)

Calculate the complex velocity in the circle plane ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit angular velocity
about its reference point.
"""
@inline wcbr(ζ,b::ConformalBody) = conj(_wcbr_conj(ζ,b))

"""
    wcv(ζ,b::ConformalBody,v::Element)

Calculate the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to a unit-strength vortex at the location of `v`,
along with its images inside the body `b`. It is presumed that the location of `v` is given in the circle plane.
"""
@inline wcv(ζ,b::ConformalBody,v::Element) = conj(_wcv_conj(ζ,b,v))

@inline _wcbt1_conj(ζ,b::ConformalBody) = induce_velocity(ζ,_set_motion_in_quiescent(b,1.0,0.0,0.0),0.0)
@inline _wcbt2_conj(ζ,b::ConformalBody) = induce_velocity(ζ,_set_motion_in_quiescent(b,0.0,1.0,0.0),0.0)
@inline _wcbr_conj(ζ,b::ConformalBody) = induce_velocity(ζ,_set_motion_in_quiescent(b,0.0,0.0,1.0),0.0)
@inline _wcv_conj(ζ,b::ConformalBody,v::Element) = induce_velocity(ζ,(_set_image_for_stationary_body(b,v,Val(true)),_unit_strength_copy(v)),0.0)

#### Derivatives of circle-plane basis velocity fields ####
"""
    dwcbt1dζ(ζ,b::ConformalBody)

Derivative with respect to ``\\zeta`` of the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline dwcbt1dζ(ζ,b::ConformalBody) = dinduce_velocity_dz(ζ,_set_motion_in_quiescent(b,1.0,0.0,0.0),0.0)

"""
    dwcbt1dζstar(ζ,b::ConformalBody)

Derivative with respect to ``\\zeta^*`` of the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline dwcbt1dζstar(ζ,b::ConformalBody) = dinduce_velocity_dzstar(ζ,_set_motion_in_quiescent(b,1.0,0.0,0.0),0.0)

"""
    dwcbt2dζ(ζ,b::ConformalBody)

Derivative with respect to ``\\zeta`` of the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline dwcbt2dζ(ζ,b::ConformalBody) = dinduce_velocity_dz(ζ,_set_motion_in_quiescent(b,0.0,1.0,0.0),0.0)

"""
    dwcbt2dζstar(ζ,b::ConformalBody)

Derivative with respect to ``\\zeta^*`` of the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline dwcbt2dζstar(ζ,b::ConformalBody) = dinduce_velocity_dzstar(ζ,_set_motion_in_quiescent(b,0.0,1.0,0.0),0.0)

"""
    dwcbrdζ(ζ,b::ConformalBody)

Derivative with respect to ``\\zeta`` of the complex velocity in the circle plane ``wc = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit angular velocity
about its reference point.
"""
@inline dwcbrdζ(ζ,b::ConformalBody) = dinduce_velocity_dz(ζ,_set_motion_in_quiescent(b,0.0,0.0,1.0),0.0)

"""
    dwcbrdζstar(ζ,b::ConformalBody)

Derivative with respect to ``\\zeta^*`` of the complex velocity in the circle plane ``wc = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit angular velocity
about its reference point.
"""
@inline dwcbrdζstar(ζ,b::ConformalBody) = dinduce_velocity_dzstar(ζ,_set_motion_in_quiescent(b,0.0,0.0,1.0),0.0)


"""
    dwcvdζ(ζ,b::ConformalBody,v::Element)

Derivative of vortex-induced basis velocity at `ζ` (in circle plane) with respect to ``\\zeta``
"""
dwcvdζ(ζ,b::Bodies.ConformalBody,v::Element) = dinduce_velocity_dz(ζ,(_set_image_for_stationary_body(b,v,Val(true)),_unit_strength_copy(v)),0.0)

"""
    dwcvdζstar(ζ,b::ConformalBody,v::Element)

Derivative of vortex-induced basis velocity at `ζ` (in circle plane) with respect to ``\\zeta^*``
"""
dwcvdζstar(ζ,b::Bodies.ConformalBody,v::Element) = dinduce_velocity_dzstar(ζ,_unit_strength_copy(v),0.0)

"""
    dwcvdζv(ζ,b::ConformalBody,v::Element)

Derivative of vortex-induced basis velocity at `ζ` (in circle plane) with respect to
the change of position of element `v`.
"""
dwcvdζv(ζ,b::Bodies.ConformalBody,v::Element) = -dinduce_velocity_dz(ζ,_unit_strength_copy(v),0.0)

"""
    dwcvdζvstar(ζ,b::ConformalBody,v::Element)

Derivative of vortex-induced basis velocity at `ζ` (in circle plane) with respect to
the change of conjugate position of element `v`.
"""
function dwcvdζvstar(ζ,b::Bodies.ConformalBody,v::Element)
  b_img = _set_image_for_stationary_body(b,v,Val(false))
  img = first(b_img.img)
  return -dinduce_velocity_dzstar(ζ,_unit_strength_copy(v),0.0) + dinduce_velocity_dz(ζ,b_img,0.0)*position(img)^2
end

#### Physical plane basis velocity fields ####


"""
    wbt1(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline wbt1(ζ,b::Bodies.ConformalBody) = conj(Bodies.transform_velocity(_wcbt1_conj(ζ,b),ζ,b))

"""
    wbt2(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to motion of body `b` with unit velocity
in the ỹ direction (of its own coordinate system).
"""
@inline wbt2(ζ,b::Bodies.ConformalBody) = conj(Bodies.transform_velocity(_wcbt2_conj(ζ,b),ζ,b))

"""
    wbr(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to motion of body `b` with unit velocity
angular velocity.
"""
@inline wbr(ζ,b::Bodies.ConformalBody) = conj(Bodies.transform_velocity(_wcbr_conj(ζ,b),ζ,b))


"""
    winf1(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to unit velocity at infinity
in the x̃ direction (of its own coordinate system).
"""
@inline winf1(ζ,b::Bodies.ConformalBody) = 1.0 .- wbt1(ζ,b)

"""
    winf2(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to unit velocity at infinity
in the ỹ direction (of its own coordinate system).
"""
@inline winf2(ζ,b::Bodies.ConformalBody) = -im .- wbt2(ζ,b)

"""
    wrinf(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to unit angular velocity at infinity
about the body's reference point.
"""
@inline wrinf(ζ,b::Bodies.ConformalBody) = conj(im*b.m.(ζ)*exp(im*b.α)) .- wbr(ζ,b)

"""
    wv(ζ,b::ConformalBody,v::Element)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to a unit-strength vortex at the location of `v`,
along with its images inside the body `b`. It is presumed that the location of `v` is given in the circle plane.
"""
@inline wv(ζ,b::Bodies.ConformalBody,v::Element) = conj(Bodies.transform_velocity(_wcv_conj(ζ,b,v),ζ,b))

#### Derivatives of physical plane basis velocity fields ####

"""
    dwvdz(ζ,b::ConformalBody,v::Element)

Derivative with respect to ``z`` of basis vortex-induced velocity ``w_v`` from vortex `v`
at physical plane mapped from location `ζ` (in circle plane).
"""
@inline dwvdz(ζ,b::Bodies.ConformalBody,v::Element) = conj(Bodies.transform_velocity(conj(_dwvdζ(ζ,b,v)),ζ,b))

"""
    dwvdzstar(ζ,b::ConformalBody,v::Element)

Derivative with respect to conjugate of ``z`` of basis vortex-induced velocity ``w_v`` from vortex `v`
at physical plane mapped from location `ζ` (in circle plane).
"""
@inline dwvdzstar(ζ,b::Bodies.ConformalBody,v::Element) = Bodies.transform_velocity(_dwvdζstar(ζ,b,v),ζ,b)

"""
    dwvdzv(ζ,b::ConformalBody,v::Element)

Derivative with respect to vortex position of basis vortex-induced velocity ``w_v`` from vortex `v`
evaluate at physical plane mapped from location `ζ` (in circle plane).
"""
@inline dwvdzv(ζ,b::Bodies.ConformalBody,v::Element) = conj(Bodies.transform_velocity(conj(_dwvdζv(ζ,b,v)),v.z,b))

"""
    dwvdzvstar(ζ,b::ConformalBody,v::Element)

Derivative with respect to conjugate vortex position of basis vortex-induced velocity ``w_v`` from vortex `v`
evaluate at physical plane mapped from location `ζ` (in circle plane).
"""
@inline dwvdzvstar(ζ,b::Bodies.ConformalBody,v::Element) = Bodies.transform_velocity(_dwvdζvstar(ζ,b,v),v.z,b)


# Derivatives of winf velocities with respect to physical plane coordinates
for w in [:inf1,:inf2,:rinf,:bt1,:bt2,:br]
  dwdz = Symbol("dw",w,"dz")
  dwdζ = Symbol("_dw",w,"dζ")
  dwdzstar = Symbol("dw",w,"dzstar")
  dwdζstar = Symbol("_dw",w,"dζstar")

  @eval function $dwdz(ζ,b::Bodies.ConformalBody)
    out = $dwdζ(ζ,b)
    return conj(Bodies.transform_velocity(conj(out),ζ,b))
  end

  @eval function $dwdzstar(ζ,b::Bodies.ConformalBody)
    out = $dwdζstar(ζ,b)
    return Bodies.transform_velocity(out,ζ,b)
  end

end

#=
function dwinf1dz(ζ,b::Bodies.ConformalBody)
  out = _dwinf1dζ(ζ,b)
  return conj(Bodies.transform_velocity(conj(out),ζ,b))
end
=#

# Derivative of wbt1, wbt2, wbr (in physical plane) with respect to zeta, zeta* (in circle plane)
for w in [:bt1,:bt2,:br]
    wc = Symbol("wc",w)
    dwcdζ = Symbol("d",wc,"dζ")
    dwdζ = Symbol("_dw",w,"dζ")
    dwcdζstar = Symbol("d",wc,"dζstar")
    dwdζstar = Symbol("_dw",w,"dζstar")

    @eval function $dwdζ(ζ,b::Bodies.ConformalBody)
        dz̃, ddz̃, _ = derivatives(ζ,b.m)
        out = Bodies.$dwcdζ(ζ,b) - ddz̃/dz̃*Bodies.$wc(ζ,b)
        return conj(Bodies.transform_velocity(conj(out),ζ,b))
    end

    @eval function $dwdζstar(ζ,b::Bodies.ConformalBody)
        out = Bodies.$dwcdζstar(ζ,b)
        return conj(Bodies.transform_velocity(conj(out),ζ,b))
    end
end

# Derivative of wv (in physical plane) with respect to zeta (in circle plane)
function _dwvdζ(ζ,b::Bodies.ConformalBody,v::Element)
  dz̃, ddz̃, _ = derivatives(ζ,b.m)
  out = dwcvdζ(ζ,b,v) - ddz̃/dz̃*wcv(ζ,b,v)
  return conj(Bodies.transform_velocity(conj(out),ζ,b))
end

# Derivative of wv (in physical plane) with respect to zeta* (in circle plane)
function _dwvdζstar(ζ,b::Bodies.ConformalBody,v::Element)
  out = dwcvdζstar(ζ,b,v)
  return conj(Bodies.transform_velocity(conj(out),ζ,b))
end

# Derivative of wv (in physical plane) with respect to vortex position (in circle plane)
function _dwvdζv(ζ,b::Bodies.ConformalBody,v::Element)
  dz̃, ddz̃, _ = derivatives(ζ,b.m)
  out = dwcvdζv(ζ,b,v)
  return conj(Bodies.transform_velocity(conj(out),ζ,b))
end

# Derivative of wv (in physical plane) with respect to conjugate vortex position (in circle plane)
function _dwvdζvstar(ζ,b::Bodies.ConformalBody,v::Element)
  dz̃, ddz̃, _ = derivatives(ζ,b.m)
  out = dwcvdζvstar(ζ,b,v)
  return conj(Bodies.transform_velocity(conj(out),ζ,b))
end

_dwinf1dζ(ζ,b::Bodies.ConformalBody) = -_dwbt1dζ(ζ,b)
_dwinf2dζ(ζ,b::Bodies.ConformalBody) = -_dwbt2dζ(ζ,b)
_dwrinfdζ(ζ,b::Bodies.ConformalBody) = -_dwbrdζ(ζ,b)

_dwinf1dζstar(ζ,b::Bodies.ConformalBody) = -_dwbt1dζstar(ζ,b)
_dwinf2dζstar(ζ,b::Bodies.ConformalBody) = -_dwbt2dζstar(ζ,b)
function _dwrinfdζstar(ζ,b::Bodies.ConformalBody)
  dz̃, ddz̃, _ = derivatives(ζ,b.m)
  return conj(im*dz̃*exp(im*b.α)) - _dwbrdζstar(ζ,b)
end

#### Some other important quantities ####

"""
    dphivdzv(ζ,b::ConformalBody,v::Element)

Calculates the derivative of the scalar potential field ``\\phi`` with respect to a change in
the body-fixed coordinates of element `v`, also accounting for the change in the image of `v`.
"""
function dphivdzv(ζ,b::Bodies.ConformalBody,v::Element)
    return conj(Bodies.transform_velocity(conj(_dphivdζv(ζ,b,v)),v.z,b))
end

function _dphivdζv(ζ,b::Bodies.ConformalBody,v::Element)
    out = conj(induce_velocity(ζ,_unit_strength_copy(v),0.0))
    b_img = _set_image_for_stationary_body(b,v,Val(false))
    img = first(b_img.img)
    out -= induce_velocity(ζ,b_img,0.0)*conj(Elements.position(img))^2
    return -0.5*out
end

"""
    wvv(targ::Element,src::Element,b::ConformalBody)

Calculate the complex velocity ``w = u - iv`` of the target element `targ` due to element `src` with unit strength,
also accounting for the images of `src`. If `targ` and `src` are the same element, then
it makes the Routh correction.
"""
@inline wvv(targ::Element,src::Element,b::Bodies.ConformalBody) = wv(targ.z,b,src) +
                                                                   conj(Bodies.transform_velocity(conj(_routh(targ.z,b,Val(targ==src))),targ.z,b))

#### Derivatives ####
"""
    d2phivdzv2(ζ,b::ConformalBody,v::Element)

Calculates the second derivative of the scalar potential field ``\\phi`` with respect to a change in
the physical-plane coordinates of element `v`, also accounting for the change in the image of `v`.
"""
function d2phivdzv2(ζ,b::ConformalBody,v::Element)
  dz̃, ddz̃, _ = derivatives(v.z,b.m)
  out = d2phivdζv2(ζ,b,v) - ddz̃/dz̃*_dphivdζv(ζ,b,v)
  out = conj(transform_velocity(conj(out),v.z,b))
  out = conj(transform_velocity(conj(out),v.z,b))
  return out
end

"""
    d2phivdzvdzvstar(ζ,b::ConformalBody,v::Element)

Calculates the second derivative of the scalar potential field ``\\phi`` with respect to a change in
the physical-plane coordinates of element `v`, also accounting for the change in the image of `v`.
"""
function d2phivdzvdzvstar(ζ,b::ConformalBody,v::Element)
  out = d2phivdζvdζvstar(ζ,b,v)
  out = conj(transform_velocity(conj(out),v.z,b))
  out = transform_velocity(out,v.z,b)
  return out
end

"""
    d2phivdζv2(ζ,b::ConformalBody,v::Element)

Calculates the second derivative of the scalar potential field ``\\phi`` with respect to a change in
the body-fixed coordinates of element `v`, also accounting for the change in the image of `v`.
"""
function d2phivdζv2(ζ,b::ConformalBody,v::Element)
    out = dinduce_velocity_dz(ζ,_unit_strength_copy(v),0.0)
    b_img = _set_image_for_stationary_body(b,v,Val(false))
    img = first(b_img.img)
    out -= 2.0*induce_velocity(ζ,b_img,0.0)*conj(position(img)^3)
    out += conj(dinduce_velocity_dz(ζ,b_img,0.0)*position(img)^4)
    return 0.5*out
end

"""
    d2phivdζvdζvstar(ζ,b::ConformalBody,v::Element)

Calculates the second derivative of the scalar potential field ``\\phi`` with respect to a change in
the body-fixed coordinates of element `v` and its conjugate, also accounting for the change in the image of `v`.
"""
d2phivdζvdζvstar(ζ,b::ConformalBody,v::Element) = 0.5*dinduce_velocity_dzstar(ζ,_unit_strength_copy(v),0.0)

"""
    dwvvdz_src(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to physical-plane position of `src`.
"""
function dwvvdz_src(targ::Element,src::Element,b::Bodies.ConformalBody)
  out = dwvvdzeta_src(targ,src,b)
  return conj(Bodies.transform_velocity(conj(out),src.z,b))
end

"""
    dwvvdz_targ(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to physical-plane position of `targ`.
"""
function dwvvdz_targ(targ::Element,src::Element,b::Bodies.ConformalBody)
  out = dwvvdzeta_targ(targ,src,b)
  return conj(Bodies.transform_velocity(conj(out),targ.z,b))
end

"""
    dwvvdzstar_src(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to the conjugate of the physical-plane position of `src`.
"""
function dwvvdzstar_src(targ::Element,src::Element,b::Bodies.ConformalBody)
  out = dwvvdzetastar_src(targ,src,b)
  return Bodies.transform_velocity(out,src.z,b)
end

"""
    dwvvdzstar_targ(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to the conjugate of the physical-plane position of `targ`.
"""
function dwvvdzstar_targ(targ::Element,src::Element,b::Bodies.ConformalBody)
  out = dwvvdzetastar_targ(targ,src,b)
  return Bodies.transform_velocity(out,targ.z,b)
end


"""
    dwvvdzeta_src(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to position of `src` in the circle plane.
"""
function dwvvdzeta_src(targ::Element,src::Element,b::Bodies.ConformalBody)
  out = dwcvdζv(targ.z,b,src)
  return conj(Bodies.transform_velocity(conj(out),targ.z,b))
end

"""
    dwvvdzetastar_src(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to conjugate position of `src` in the circle plane.
"""
function dwvvdzetastar_src(targ::Element,src::Element,b::Bodies.ConformalBody)
  out = dwcvdζvstar(targ.z,b,src)
  return conj(Bodies.transform_velocity(conj(out),targ.z,b))
end

"""
    dwvvdzeta_targ(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to position of `targ` in the circle plane.
"""
dwvvdzeta_targ(targ::Element,src::Element,b::Bodies.ConformalBody) =
                            _dwvdζ(targ.z,b,src) +
                            conj(Bodies.transform_velocity(conj(_drouthdzeta(targ.z,b,Val(targ==src))),targ.z,b))

"""
    dwvvdzetastar_targ(targ::Element,src::Element,b::ConformalBody)

Calculate the derivative of the complex velocity ``w = u - iv`` of
induced on the target element `targ` by the element `src` with unit strength,
with respect to conjugate position of `targ` in the circle plane.
"""
dwvvdzetastar_targ(targ::Element,src::Element,b::Bodies.ConformalBody) =
                    _dwvvdzetastar_targ(targ,src,b,Val(targ==src))


_dwvvdzetastar_targ(targ,src,b,::Val{true}) = zero(targ.z)

function _dwvvdzetastar_targ(targ,src,b,::Val{false})
  out = dwcvdζstar(targ.z,b,src)
  return conj(Bodies.transform_velocity(conj(out),targ.z,b))
end

#### Routh correction and its derivatives ####

function _routh(ζ,b,::Val{true})
    dz̃, ddz̃, _ = derivatives(ζ,b.m)
    return -ddz̃/(4π*im*dz̃)
end

_routh(ζ,b,::Val{false}) = zero(ζ)

function _drouthdzeta(ζ,b,::Val{true})
    dz̃, ddz̃, dddz̃ = derivatives(ζ,b.m)
    return -(dddz̃ - 2ddz̃^2/dz̃)/(4π*im*dz̃)
end

_drouthdzeta(ζ,b,::Val{false}) = zero(ζ)



#### Some utilities ####


function _set_motion_in_quiescent(b_orig,Ũ,Ṽ,Ω)
    b = deepcopy(b_orig)
    clear_images!(b)
    motion = RigidBodyMotion((Ũ+im*Ṽ)*exp(im*b.α), Ω)
    Bodies.enforce_no_flow_through!(b, motion, () , 0.0)
    return b
end

function _set_image_for_stationary_body(b_orig,v::T,atorigin) where {T<:Element}
    b = deepcopy(b_orig)
    clear_images!(b)
    motion = RigidBodyMotion(0.0,0.0)
    sys = _preimage_system(v,atorigin)
    Bodies.enforce_no_flow_through!(b, motion, sys , 0.0)
    return b
end

_preimage_system(v,::Val{true}) = (_unit_strength_copy(v),_unit_at_origin(v)) # put an image at origin, too
_preimage_system(v,::Val{false}) = (_unit_strength_copy(v),) # for no image at origin
_unit_strength_copy(v::Blob{T}) where {T} = Blob{T}(v.z,1.0,v.δ)
_unit_strength_copy(v::Point{T}) where {T} = Point{T}(v.z,1.0)
_unit_at_origin(::Union{Blob{T},Point{T}}) where {T} = Point{T}(complex(Inf),-1.0)
