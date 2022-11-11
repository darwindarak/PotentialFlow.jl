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
    force(v::Vector{Element},b::ConformalBody) -> Tuple{Float64,Float64,Float64}

Return the force and moment due to the vortex elements in `v` and any motion in `b`.
Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function force(v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}

        mr, fx, fy = 0.0, 0.0, 0.0

        Ma = addedmass(b)

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇
        Uvec = [Ω,Ũ,Ṽ]

        Ũ̇, Ṽ̇ = reim(b.c̈*exp(-im*b.α))
        Ω̇ = b.α̈
        U̇vec = [Ω̇,Ũ̇,Ṽ̇]

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

            fx += Γj*FUv_x(vj,b)*Uvec
            fy += Γj*FUv_y(vj,b)*Uvec
            mr += Γj*FUv_r(vj,b)*Uvec
        end

        P̃U = Ma*Uvec
        fx += Ω*P̃U[3]
        fy -= Ω*P̃U[2]
        mr += Ṽ*P̃U[2]-Ũ*P̃U[3]

        f = -Ma*U̇vec
        mr += f[1]
        fx += f[2]
        fy += f[3]

        return fx, fy, mr
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


"""
    dfdζv(l::Integer,v::Vector{Element},b::ConformalBody) -> Tuple{ComplexF64,ComplexF64,ComplexF64}

Return the derivatives with respect to the change of circle-plane ``\\zeta`` position of the vortex with index `l`
of the force and moment,
due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dfdζv(l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}

  dfx, dfy, dmr = dfdzv(l,v,b;kwargs...)
  dz̃, ddz̃, _ = derivatives(v[l].z,b.m)
  dfx *= dz̃*exp(im*b.α)
  dfy *= dz̃*exp(im*b.α)
  dmr *= dz̃*exp(im*b.α)

  return dfx, dfy, dmr
end

"""
    dfdzv(l::Integer,v::Vector{Element},b::ConformalBody) -> Tuple{ComplexF64,ComplexF64,ComplexF64}

Return the derivatives with respect to the change of physical-plane ``z`` position of the vortex with index `l`
of the force and moment,
due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dfdzv(l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
        @assert l <= length(v) "Index exceeds length of vector of elements"

        dfx, dfy, dmr = complex(0.0), complex(0.0), complex(0.0)

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇
        Uvec = [Ω,Ũ,Ṽ]

        vl = v[l]
        Γl = circulation(vl)

        for (j,vj) in enumerate(v)
            Γj  = circulation(vj)
            dfx -= Γj*Γl*dFvv_xdzv(vl,vj,b;kwargs...)
            dfy -= Γj*Γl*dFvv_ydzv(vl,vj,b;kwargs...)
            dmr -= Γj*Γl*dFvv_rdzv(vl,vj,b;kwargs...)

        end
        dfx += Γl*dFUv_xdzv(vl,b;kwargs...)*Uvec
        dfy += Γl*dFUv_ydzv(vl,b;kwargs...)*Uvec
        dmr += Γl*dFUv_rdzv(vl,b;kwargs...)*Uvec

        return dfx,dfy,dmr
end

"""
    dfdΓv(l::Integer,v::Vector{Element},b::ConformalBody) -> Tuple{Float64,Float64,Float64}

Return the derivatives of the force and moment with respect to
the change of circulation of the vortex with index `l`, due to the vortex elements in `v` and any motion in `b`. Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dfdΓv(l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}

        dmr, dfx, dfy = 0.0, 0.0, 0.0

        Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
        Ω = b.α̇
        Uvec = [Ω,Ũ,Ṽ]

        vl = v[l]
        Γl = circulation(vl)

        for (j,vj) in enumerate(v)
            Γj  = circulation(vj)
            dfx -= Γj*Fvv_x(vl,vj,b;kwargs...)
            dfy -= Γj*Fvv_y(vl,vj,b;kwargs...)
            dmr -= Γj*Fvv_r(vl,vj,b;kwargs...)
        end
        dfx += FUv_x(vl,b)*Uvec
        dfy += FUv_y(vl,b)*Uvec
        dmr += FUv_r(vl,b)*Uvec

        return dfx, dfy, dmr
end

"""
    dfdU(ζ,l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative of the force and moment
due to the vortex elements in `v` and any motion in `b`, with respect to
the change of rigid-body motion component with index `l`. Note that
these components are index as follows: ``[\\Omega,\\tilde{U}_r,\\tilde{V}_r]``.
Note that it is presumed that the elements `v` are given in
the circle plane.
"""
function dfdU(l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}

  dmr, dfx, dfy = 0.0, 0.0, 0.0

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇

  for (j,vj) in enumerate(v)
      Γj  = circulation(vj)
      dfx += Γj*FUv_x(vj,b)[l]
      dfy += Γj*FUv_y(vj,b)[l]
      dmr += Γj*FUv_r(vj,b)[l]
  end
  dfxU, dfyU, dmrU = _dfdU(b,Val(l);kwargs...)

  return dfx+dfxU, dfy+dfyU, dmr+dmrU
end

# derivative wrt Ω
function _dfdU(b,::Val{1};kwargs...)
  dmr, dfx, dfy = 0.0, 0.0, 0.0

  Ma = addedmass(b)

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇
  Uvec = [Ω,Ũ,Ṽ]
  P̃U = Ma*Uvec

  dfx =  P̃U[3] + Ω*Ma[3,1]
  dfy = -P̃U[2] - Ω*Ma[2,1]
  dmr = Ṽ*Ma[2,1] - Ũ*Ma[3,1]

  return dfx, dfy, dmr
end

# derivative wrt Ũ
function _dfdU(b,::Val{2};kwargs...)
  dmr, dfx, dfy = 0.0, 0.0, 0.0

  Ma = addedmass(b)

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇
  Uvec = [Ω,Ũ,Ṽ]
  P̃U = Ma*Uvec

  dfx =  Ω*Ma[3,2]
  dfy = -Ω*Ma[2,2]
  dmr = -P̃U[3] + Ṽ*Ma[2,2] - Ũ*Ma[3,2]

  return dfx, dfy, dmr
end

# derivative wrt Ṽ
function _dfdU(b,::Val{3};kwargs...)
  dmr, dfx, dfy = 0.0, 0.0, 0.0

  Ma = addedmass(b)

  Ũ, Ṽ = reim(b.ċ*exp(-im*b.α))
  Ω = b.α̇
  Uvec = [Ω,Ũ,Ṽ]
  P̃U = Ma*Uvec

  dfx =  Ω*Ma[3,3]
  dfy = -Ω*Ma[2,3]
  dmr = P̃U[2] + Ṽ*Ma[2,3] - Ũ*Ma[3,3]

  return dfx, dfy, dmr
end

"""
    dfdUdot(l::Integer,v::Vector{Element},b::ConformalBody)

Return the derivative of the force and moment
due to the vortex elements in `v` and any motion in `b`, with respect to
the change of rigid-body acceleration component with index `l`. Note that
these components are index as follows: ``[\\dot{\\Omega},\\dot{\\tilde{U}}_r,\\dot{\\tilde{V}}_r]``
"""
function dfdUdot(l::Integer,v::Vector{T},b::Bodies.ConformalBody;kwargs...) where {T<:Element}
  Ma = addedmass(b)
  return -Ma[2,l], -Ma[3,l], -Ma[1,l]
end


#=
function P(ζ,v::Element,b::Bodies.ConformalBody)
    return abs2.(wv(ζ,b,v)) + 4.0*real(dphivdzv(ζ,b,v)*conj(wvv(v,v,b)))
end
=#

function Πvv(ζ,targ::Element,src::Element,b::Bodies.ConformalBody;kwargs...)
    return real(wv(ζ,b,targ;kwargs...).*conj(wv(ζ,b,src;kwargs...))) +
           2.0*real(dphivdzv(ζ,b,targ;kwargs...)*conj(wvv(targ,src,b;kwargs...))) +
           2.0*real(dphivdzv(ζ,b,src)*conj(wvv(src,targ,b;kwargs...)))
end

dΠvvdzv(ζ,targ::Element,src::Element,b::Bodies.ConformalBody;kwargs...) = _dΠvvdzv_targ(ζ,targ,src,b;kwargs...)

# This computes the derivative of Πvv with respect to the
# target's position. Note that the derivative with respect to
# the source's position is simply the same with targ and src swapped.
function _dΠvvdzv_targ(ζ,targ::Element,src::Element,b::Bodies.ConformalBody;kwargs...)
    out = dwvdzv(ζ,b,targ).*conj(wv(ζ,b,src;kwargs...)) +
          2.0*d2phivdzv2(ζ,b,targ)*conj(wvv(targ,src,b;kwargs...)) +
          2.0*dphivdzv(ζ,b,targ)*conj(dwvvdzstar_targ(targ,src,b)) +
          2.0*dphivdzv(ζ,b,src)*conj(dwvvdzstar_src(src,targ,b))
    out += conj(dwvdzvstar(ζ,b,targ).*conj(wv(ζ,b,src;kwargs...)) +
          2.0*d2phivdzvdzvstar(ζ,b,targ)*conj(wvv(targ,src,b;kwargs...)) +
          2.0*dphivdzv(ζ,b,targ)*conj(dwvvdz_targ(targ,src,b;kwargs...)) +
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

for (wtype,c,comp) in [(:inf1,:U,:x),(:inf2,:V,:y),(:rinf,:Ω,:r)]
  Πfcn = Symbol("Π",c,"v")
  dΠfcn = Symbol("dΠ",c,"vdzv")

  Ffcn = Symbol("Fvv_",comp)
  dFfcn = Symbol("dFvv_",comp,"dzv")

  winf = Symbol("w",wtype)
  dwinfdz = Symbol("d",winf,"dz")
  dwinfdzstar = Symbol("d",winf,"dzstar")

  # ΠUv, ΠVv, ΠΩv
  # switched from 4 to 2 in front of second term
  @eval function $Πfcn(ζ,v::Element,b::Bodies.ConformalBody;kwargs...)
      return -real(conj($winf(ζ,b)).*wv(ζ,b,v;kwargs...)) -
             2.0*real(dphivdzv(ζ,b,v)*conj($winf(v.z,b)))
  end

  # dΠUvdz, dΠVvdz, dΠΩvdz
  # switched from 4 to 2
  @eval function $dΠfcn(ζ,v::Element,b::Bodies.ConformalBody)
    winf_ζ_star = conj($winf(ζ,b))
    winf_v_star = conj($winf(v.z,b))
    dwinfstar_v = conj($dwinfdzstar(v.z,b))
    dwinfstar_vstar = conj($dwinfdz(v.z,b))
    out = winf_ζ_star.*dwvdzv(ζ,b,v)+2.0*winf_v_star*d2phivdzv2(ζ,b,v)+2.0*dphivdzv(ζ,b,v)*dwinfstar_v
    out += conj(winf_ζ_star.*dwvdzvstar(ζ,b,v)+2.0*winf_v_star*d2phivdzvdzvstar(ζ,b,v)+2.0*dphivdzv(ζ,b,v)*dwinfstar_vstar)
    return -0.5*out
    return -out
  end

  # Fvv_x, Fvv_y, Fvv_r
  @eval function $Ffcn(targ::Element,src::Element,b::ConformalBody;kwargs...)
    winf_targ = $winf(targ.z,b)
    winf_src = $winf(src.z,b)
    return imag(winf_targ*conj(wvv(targ,src,b;kwargs...))) + imag(winf_src*conj(wvv(src,targ,b;kwargs...)))
  end

  # dFvv_xdzv, dFvv_ydzv, dFvv_rdzv (with respect to target)
  @eval function $dFfcn(targ::Element,src::Element,b::ConformalBody;kwargs...)
    winf_targ = $winf(targ.z,b)
    dwinf_targ = $dwinfdz(targ.z,b)
    dwinfstar_targ = $dwinfdzstar(targ.z,b)
    wvv_targ = wvv(targ,src,b;kwargs...)
    winf_src = $winf(src.z,b)

    out = dwinf_targ*conj(wvv_targ)+ winf_targ*conj(dwvvdzstar_targ(targ,src,b))
    out -= conj(dwinfstar_targ)*wvv_targ + conj(winf_targ)*dwvvdz_targ(targ,src,b;kwargs...)
    out += winf_src*conj(dwvvdzstar_src(src,targ,b))
    out -= conj(winf_src)*dwvvdz_src(src,targ,b)
    return -0.5im*out
  end
end


# Note that these rely on unit_impulse of the vortex, which only
# strictly computes the impulse for the image of equal/opposite strength.

function FUv_x(v::Element,b::ConformalBody)
  Pv = unit_impulse(v,b)
  ζv = Elements.position(v)
  return transpose([imag(Pv)+imag(winf1(ζv,b)*conj(wrinf(ζv,b))),
                    0.0,
                    imag(winf1(ζv,b)*conj(winf2(ζv,b)))])
end

function FUv_y(v::Element,b::ConformalBody)
  Pv = unit_impulse(v,b)
  ζv = Elements.position(v)
  return transpose([-real(Pv)+imag(winf2(ζv,b)*conj(wrinf(ζv,b))),
                     imag(winf2(ζv,b)*conj(winf1(ζv,b))),
                     0.0])
end

function FUv_r(v::Element,b::ConformalBody)
  Pv = unit_impulse(v,b)
  ζv = Elements.position(v)
  return transpose([0.0,
                    -imag(Pv)+imag(wrinf(ζv,b)*conj(winf1(ζv,b))),
                     real(Pv)+imag(wrinf(ζv,b)*conj(winf2(ζv,b)))])
end

function dFUv_xdzv(v::Element,b::ConformalBody)
  ζv = Elements.position(v)
  winf_v = winf1(ζv,b)
  winf_v_star = conj(winf_v)
  return transpose([dPvydzv(v,b)-0.5im*(dwinf1dz(ζv,b)*conj(wrinf(ζv,b))-winf_v_star*dwrinfdz(ζv,b)-conj(dwinf1dzstar(ζv,b))*wrinf(ζv,b)+winf_v*conj(dwrinfdzstar(ζv,b))),
                    0.0,
                    -0.5im*(dwinf1dz(ζv,b)*conj(winf2(ζv,b))-winf_v_star*dwinf2dz(ζv,b))])
end

function dFUv_ydzv(v::Element,b::ConformalBody)
  ζv = Elements.position(v)
  winf_v = winf2(ζv,b)
  winf_v_star = conj(winf_v)
  return transpose([-dPvxdzv(v,b) - 0.5im*(dwinf2dz(ζv,b)*conj(wrinf(ζv,b))-winf_v_star*dwrinfdz(ζv,b)-conj(dwinf2dzstar(ζv,b))*wrinf(ζv,b)+winf2(ζv,b)*conj(dwrinfdzstar(ζv,b))),
                    -0.5im*(dwinf2dz(ζv,b)*conj(winf1(ζv,b))-winf_v_star*dwinf1dz(ζv,b)),
                    0.0])
end

function dFUv_rdzv(v::Element,b::ConformalBody)
  ζv = Elements.position(v)
  return transpose([0.0,
                    -dPvydzv(v,b)-0.5im*(dwrinfdz(ζv,b)*conj(winf1(ζv,b))+wrinf(ζv,b)*conj(dwinf1dzstar(ζv,b))-conj(dwrinfdzstar(ζv,b))*winf1(ζv,b)-conj(wrinf(ζv,b))*dwinf1dz(ζv,b)),
                     dPvxdzv(v,b)-0.5im*(dwrinfdz(ζv,b)*conj(winf2(ζv,b))+wrinf(ζv,b)*conj(dwinf2dzstar(ζv,b))-conj(dwrinfdzstar(ζv,b))*winf2(ζv,b)-conj(wrinf(ζv,b))*dwinf2dz(ζv,b))])
end
