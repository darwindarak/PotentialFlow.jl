#=
 Routines that provide the unit basis fields for a given body
 Note that the unit velocity fields are provided in standard
 complex form, w = u - iv
=#

export Fbt1, Fbt2, Fbr, Fv, wbt1, wbt2, wbr, winf1, winf2, wrinf, wv,
        ŵbt1, ŵbt2, ŵbr, ŵv

"""
    Fbt1(ζ,b::ConformalBody)

Calculate the complex potential at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline Fbt1(ζ,b::Bodies.ConformalBody) = Elements.complexpotential(ζ,_set_motion_in_quiescent(b,1.0,0.0,0.0))

"""
    Fbt2(ζ,b::ConformalBody)

Calculate the complex potential at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the ỹ direction (of its own coordinate system).
"""
@inline Fbt2(ζ,b::Bodies.ConformalBody) = Elements.complexpotential(ζ,_set_motion_in_quiescent(b,0.0,1.0,0.0))

"""
    Fbr(ζ,b::ConformalBody)

Calculate the complex potential at `ζ` in the circle plane due to motion of body `b` with unit angular velocity
about its reference point
"""
@inline Fbr(ζ,b::Bodies.ConformalBody) = Elements.complexpotential(ζ,_set_motion_in_quiescent(b,0.0,0.0,1.0))

"""
    Fv(ζ,b::ConformalBody,v::Element)

Calculate the complex potential at `ζ` in the circle plane due to a unit-strength vortex at the location of `v`,
along with its images inside the body `b`. It is presumed that the location of `v` is given in the circle plane.
"""
@inline Fv(ζ,b::Bodies.ConformalBody,v::Element) = Elements.complexpotential(ζ,(_set_image_for_stationary_body(b,v,Val(true)),_unit_strength_copy(v)))

"""
    ŵbt1(ζ,b::ConformalBody)

Calculate the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline ŵbt1(ζ,b::Bodies.ConformalBody) = conj(_ŵbt1_conj(ζ,b))

"""
    ŵbt2(ζ,b::ConformalBody)

Calculate the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline ŵbt2(ζ,b::Bodies.ConformalBody) = conj(_ŵbt2_conj(ζ,b))

"""
    ŵbr(ζ,b::ConformalBody)

Calculate the complex velocity in the circle plane ``ŵ = û-iv̂`` at `ζ` in the circle plane due to motion of body `b` with unit angular velocity
about its reference point.
"""
@inline ŵbr(ζ,b::Bodies.ConformalBody) = conj(_ŵbr_conj(ζ,b))

"""
    ŵv(ζ,b::ConformalBody,v::Element)

Calculate the complex velocity ``ŵ = û-iv̂`` at `ζ` in the circle plane due to a unit-strength vortex at the location of `v`,
along with its images inside the body `b`. It is presumed that the location of `v` is given in the circle plane.
"""
@inline ŵv(ζ,b::Bodies.ConformalBody,v::Element) = conj(_ŵv_conj(ζ,b,v))

@inline _ŵbt1_conj(ζ,b::Bodies.ConformalBody) = induce_velocity(ζ,_set_motion_in_quiescent(b,1.0,0.0,0.0),0.0)
@inline _ŵbt2_conj(ζ,b::Bodies.ConformalBody) = induce_velocity(ζ,_set_motion_in_quiescent(b,0.0,1.0,0.0),0.0)
@inline _ŵbr_conj(ζ,b::Bodies.ConformalBody) = induce_velocity(ζ,_set_motion_in_quiescent(b,0.0,0.0,1.0),0.0)
@inline _ŵv_conj(ζ,b::Bodies.ConformalBody,v::Element) = induce_velocity(ζ,(_set_image_for_stationary_body(b,v,Val(true)),_unit_strength_copy(v)),0.0)


"""
    wbt1(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to motion of body `b` with unit velocity
in the x̃ direction (of its own coordinate system).
"""
@inline wbt1(ζ,b::Bodies.ConformalBody) = conj(Bodies.transform_velocity(_ŵbt1_conj(ζ,b),ζ,b))

"""
    wbt2(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to motion of body `b` with unit velocity
in the ỹ direction (of its own coordinate system).
"""
@inline wbt2(ζ,b::Bodies.ConformalBody) = conj(Bodies.transform_velocity(_ŵbt2_conj(ζ,b),ζ,b))

"""
    wbr(ζ,b::ConformalBody)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to motion of body `b` with unit velocity
angular velocity.
"""
@inline wbr(ζ,b::Bodies.ConformalBody) = conj(Bodies.transform_velocity(_ŵbr_conj(ζ,b),ζ,b))


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
@inline wrinf(ζ,b::Bodies.ConformalBody) = -im*conj(b.m.(ζ)) .- wbr(ζ,b)

"""
    wv(ζ,b::ConformalBody,v::Element)

Calculate the complex velocity ``w = u-iv`` in the physical plane at
the location mapped from `ζ` in the circle plane due to a unit-strength vortex at the location of `v`,
along with its images inside the body `b`. It is presumed that the location of `v` is given in the circle plane.
"""
@inline wv(ζ,b::Bodies.ConformalBody,v::Element) = conj(Bodies.transform_velocity(_ŵv_conj(ζ,b,v),ζ,b))



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
