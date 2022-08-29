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

        out -= Ũ̇*ΦU(ζ,b)
        out -= Ṽ̇*ΦV(ζ,b)
        out -= Ω̇*ΦΩ(ζ,b)


        return out
end


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

# This computes the derivative of Πvv with respect to the
# target's position. Note that the derivative with respect to
# the source's position is simply the same with targ and src swapped.
function dΠvvdzv(ζ,targ::Element,src::Element,b::Bodies.ConformalBody)
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
    out = winf_ζ_star*dwvdzv(ζ,b,v)+4.0*winf_v_star*d2phivdzv2(ζ,b,v)+4.0*dphivdzv(ζ,b,v)*dwinfstar_v
    out += conj(winf_ζ_star*dwvdzvstar(ζ,b,v)+4.0*winf_v_star*d2phivdzvdzvstar(ζ,b,v)+4.0*dphivdzv(ζ,b,v)*dwinfstar_vstar)
    return -0.5*out
    return -out
  end
end
