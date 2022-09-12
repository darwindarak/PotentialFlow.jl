"""
    pressure(ζ,v::Vector{Element},b::ConformalBody)

Return the pressure at `ζ` in the circle plane (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`.
Note that it is presumed that the elements `v` are given in
the circle plane.
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

        out -= Ũ̇*ΦU(ζ,b)
        out -= Ṽ̇*ΦV(ζ,b)
        out -= Ω̇*ΦΩ(ζ,b)


        return out
end

"""
    dpdxv(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative with respect to the change of ``x`` position of the vortex with index `l`
of the pressure at `ζ` in the circle plane (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dpdxv(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
  dp = dpdzv(ζ,l,v,b;kwargs...)
  return 2real(dp)
end

"""
    dpdyv(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative with respect to the change of ``y`` position of the vortex with index `l`
of the pressure at `ζ` in the circle plane (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dpdyv(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
  dp = dpdzv(ζ,l,v,b;kwargs...)
  return -2imag(dp)
end

"""
    dpdζv(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative with respect to the change of circle-plane ``\\zeta`` position of the vortex with index `l`
of the pressure at `ζ` in the circle plane (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dpdζv(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}

  out = dpdzv(ζ,l,v,b;kwargs...)
  dz̃, ddz̃, _ = derivatives(v[l].z,b.m)
  out *= dz̃*exp(im*b.α)
  return out
end

"""
    dpdzv(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative with respect to the change of physical-plane ``z`` position of the vortex with index `l`
of the pressure at `ζ` in the circle plane (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dpdzv(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
        @assert l <= length(v) "Index exceeds length of vector of elements"

        out = zero(ζ)

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇

        vl = v[l]
        Γl = circulation(vl)

        for (j,vj) in enumerate(v)
            Γj  = circulation(vj)
            out -= Γj*Γl*dΠvvdzv(ζ,vl,vj,b;kwargs...)
        end
        out -= Γl*Ũ*dΠUvdzv(ζ,vl,b;kwargs...)
        out -= Γl*Ṽ*dΠVvdzv(ζ,vl,b;kwargs...)
        out -= Γl*Ω*dΠΩvdzv(ζ,vl,b;kwargs...)

        return out
end

"""
    dpdΓv(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative of the pressure at `ζ` (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`, with respect to
the change of circulation of the vortex with index `l`.  Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dpdΓv(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
        @assert l <= length(v) "Index exceeds length of vector of elements"

        out = real(zero(ζ))

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇

        vl = v[l]

        for (j,vj) in enumerate(v)
            Γj  = circulation(vj)
            out -= Γj*Πvv(ζ,vl,vj,b;kwargs...)
        end
        out -= Ũ*ΠUv(ζ,vl,b;kwargs...)
        out -= Ṽ*ΠVv(ζ,vl,b;kwargs...)
        out -= Ω*ΠΩv(ζ,vl,b;kwargs...)

        return out
end

"""
    dpdU(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative of the pressure at `ζ` (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`, with respect to
the change of rigid-body motion component with index `l`. Note that
these components are index as follows: ``[\\Omega,\\tilde{U}_r,\\tilde{V}_r]``.
Note that it is presumed that the elements `v` are given in
the circle plane.
"""
dpdU(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element} =
    _dpdU(ζ,v,b,Val(l);kwargs...)


function _dpdU(ζ,v,b,::Val{1};kwargs...)
  out = real(zero(ζ))

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇

  for (j,vj) in enumerate(v)
      Γj  = circulation(vj)
      out -= Γj*ΠΩv(ζ,vj,b;kwargs...)
  end
  out -= Ω*ΠΩΩ(ζ,b)
  out -= Ũ*ΠUΩ(ζ,b)
  out -= Ṽ*ΠVΩ(ζ,b)

  return out
end

function _dpdU(ζ,v,b,::Val{2};kwargs...)
  out = real(zero(ζ))

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇

  for (j,vj) in enumerate(v)
      Γj  = circulation(vj)
      out -= Γj*ΠUv(ζ,vj,b;kwargs...)
  end
  out -= Ω*ΠUΩ(ζ,b)
  out -= Ũ*ΠUU(ζ,b)
  out -= Ṽ*ΠUV(ζ,b)

  return out
end

function _dpdU(ζ,v,b,::Val{3};kwargs...)
  out = real(zero(ζ))

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇

  for (j,vj) in enumerate(v)
      Γj  = circulation(vj)
      out -= Γj*ΠVv(ζ,vj,b;kwargs...)
  end
  out -= Ω*ΠVΩ(ζ,b)
  out -= Ũ*ΠUV(ζ,b)
  out -= Ṽ*ΠVV(ζ,b)

  return out
end

"""
    dpdUdot(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative of the pressure at `ζ` (which can be an array of points),
due to the vortex elements in `v` and any motion in `b`, with respect to
the change of rigid-body acceleration component with index `l`. Note that
these components are index as follows: ``[\\dot{\\Omega},\\dot{\\tilde{U}}_r,\\dot{\\tilde{V}}_r]``
"""
dpdUdot(ζ,l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element} =
    _dpdUdot(ζ,v,b,Val(l);kwargs...)


_dpdUdot(ζ,v,b,::Val{1};kwargs...) = -ΦΩ(ζ,b)
_dpdUdot(ζ,v,b,::Val{2};kwargs...) = -ΦU(ζ,b)
_dpdUdot(ζ,v,b,::Val{3};kwargs...) = -ΦV(ζ,b)


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

dΠvvdzv(ζ,targ::Element,src::Element,b::Bodies.ConformalBody) = _dΠvvdzv_targ(ζ,targ,src,b)

# This computes the derivative of Πvv with respect to the
# target's position. Note that the derivative with respect to
# the source's position is simply the same with targ and src swapped.
function _dΠvvdzv_targ(ζ,targ::Element,src::Element,b::Bodies.ConformalBody)
    out = dwvdzv(ζ,b,targ).*conj(wv(ζ,b,src)) +
          2.0*d2phivdzv2(ζ,b,targ)*conj(wvv(targ,src,b)) +
          2.0*dphivdzv(ζ,b,targ)*conj(dwvvdzstar_targ(targ,src,b)) +
          2.0*dphivdzv(ζ,b,src)*conj(dwvvdzstar_src(src,targ,b))
    out += conj(dwvdzvstar(ζ,b,targ).*conj(wv(ζ,b,src)) +
          2.0*d2phivdzvdzvstar(ζ,b,targ)*conj(wvv(targ,src,b)) +
          2.0*dphivdzv(ζ,b,targ)*conj(dwvvdz_targ(targ,src,b)) +
          2.0*dphivdzv(ζ,b,src)*conj(dwvvdz_src(src,targ,b)))
    return 0.5*out
end

ΠUU(ζ,b::Bodies.ConformalBody) = abs2.(winf1(ζ,b)) .- 1.0
ΠVV(ζ,b::Bodies.ConformalBody) = abs2.(winf2(ζ,b)) .- 1.0
ΠΩΩ(ζ,b::Bodies.ConformalBody) = abs2.(wrinf(ζ,b)) .- abs2.(b.m.(ζ))

ΠUV(ζ,b::Bodies.ConformalBody) = real(winf1(ζ,b).*conj(winf2(ζ,b)))
ΠUΩ(ζ,b::Bodies.ConformalBody) = real(winf1(ζ,b).*conj(wrinf(ζ,b))) + imag(b.m.(ζ))
ΠVΩ(ζ,b::Bodies.ConformalBody) = real(winf2(ζ,b).*conj(wrinf(ζ,b))) - real(b.m.(ζ))

ΦU(ζ,b::Bodies.ConformalBody) = real(Fbt1(ζ,b))
ΦV(ζ,b::Bodies.ConformalBody) = real(Fbt2(ζ,b))
ΦΩ(ζ,b::Bodies.ConformalBody) = real(Fbr(ζ,b))

# Coupled body motion - vortex terms

for (wtype,c) in [(:inf1,:U),(:inf2,:V),(:rinf,:Ω)]
  Πfcn = Symbol("Π",c,"v")
  dΠfcn = Symbol("dΠ",c,"vdzv")
  winf = Symbol("w",wtype)
  dwinfdz = Symbol("d",winf,"dz")
  dwinfdzstar = Symbol("d",winf,"dzstar")

  # ΠUv, ΠVv, ΠΩv
  @eval function $Πfcn(ζ,v::Element,b::Bodies.ConformalBody)
      return -real(conj($winf(ζ,b)).*wv(ζ,b,v)) -
             4.0*real(dphivdzv(ζ,b,v)*conj($winf(v.z,b)))
  end

  # dΠUvdz, dΠVvdz, dΠΩvdz
  @eval function $dΠfcn(ζ,v::Element,b::Bodies.ConformalBody)
    winf_ζ_star = conj($winf(ζ,b))
    winf_v_star = conj($winf(v.z,b))
    dwinfstar_v = conj($dwinfdzstar(v.z,b))
    dwinfstar_vstar = conj($dwinfdz(v.z,b))
    out = winf_ζ_star.*dwvdzv(ζ,b,v)+4.0*winf_v_star*d2phivdzv2(ζ,b,v)+4.0*dphivdzv(ζ,b,v)*dwinfstar_v
    out += conj(winf_ζ_star.*dwvdzvstar(ζ,b,v)+4.0*winf_v_star*d2phivdzvdzvstar(ζ,b,v)+4.0*dphivdzv(ζ,b,v)*dwinfstar_vstar)
    return -0.5*out
    return -out
  end
end

#### Force and moment ####
function force(v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
        fx = 0.0
        fy = 0.0
        mr = 0.0

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇
        Ũ̇, Ṽ̇ = reim(b.c̈*exp(-im*b.α))
        Ω̇ = b.α̈

        for (j,vj) in enumerate(v)
            zj,Γj  = Elements.position(vj), circulation(vj)
            fx -= 0.5*Γj^2*Fvv_x(vj,vj,b;kwargs...)
            fy -= 0.5*Γj^2*Fvv_y(vj,vj,b;kwargs...)
            mr -= 0.5*Γj^2*Fvv_r(vj,vj,b;kwargs...)
            for vk in v[1:j-1]
                zk,Γk  = Elements.position(vk), circulation(vk)
                fx -= Γj*Γk*Fvv_x(vj,vk,b;kwargs...)
                fy -= Γj*Γk*Fvv_y(vj,vk,b;kwargs...)
                mr -= Γj*Γk*Fvv_r(vj,vk,b;kwargs...)
            end

            #out -= Γj*Ũ*ΠUv(ζ,vj,b;kwargs...)
            #out -= Γj*Ṽ*ΠVv(ζ,vj,b;kwargs...)
            #out -= Γj*Ω*ΠΩv(ζ,vj,b;kwargs...)
        end
        #=
        out -= 0.5*Ũ^2*ΠUU(ζ,b)
        out -= 0.5*Ṽ^2*ΠVV(ζ,b)
        out -= 0.5*Ω^2*ΠΩΩ(ζ,b)
        out -= Ũ*Ṽ*ΠUV(ζ,b)
        out -= Ũ*Ω*ΠUΩ(ζ,b)
        out -= Ṽ*Ω*ΠVΩ(ζ,b)

        out -= Ũ̇*ΦU(ζ,b)
        out -= Ṽ̇*ΦV(ζ,b)
        out -= Ω̇*ΦΩ(ζ,b)
        =#

        return fx, fy, mr
end

# Fvv
function Fvv_x(targ::Element,src::Element,b::Bodies.ConformalBody)
  ζtarg = Elements.position(targ)
  ζsrc = Elements.position(src)
  return imag(winf1(ζtarg,b)*conj(wvv(targ,src,b))) + imag(winf1(ζsrc,b)*conj(wvv(src,targ,b)))
end

function Fvv_y(targ::Element,src::Element,b::Bodies.ConformalBody)
  ζtarg = Elements.position(targ)
  ζsrc = Elements.position(src)
  return imag(winf2(ζtarg,b)*conj(wvv(targ,src,b))) + imag(winf2(ζsrc,b)*conj(wvv(src,targ,b)))
end

function Fvv_r(targ::Element,src::Element,b::Bodies.ConformalBody)
  ζtarg = Elements.position(targ)
  ζsrc = Elements.position(src)
  return imag(wrinf(ζtarg,b)*conj(wvv(targ,src,b))) + imag(wrinf(ζsrc,b)*conj(wvv(src,targ,b)))
end
