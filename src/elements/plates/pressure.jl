import ..Properties: @property
import ..Vortex

@property begin
    signature = induce_acc(targ::Target, targvel::Target, src::Source, srcvel::Source)
    preallocator = allocate_acc
    stype = ComplexF64
end

function induce_acc(z::ComplexF64, ż::ComplexF64,
                    p::Vortex.Point, point_vel::ComplexF64)
    -circulation(p)*Points.cauchy_kernel((z - p.z)^2)*conj(ż - point_vel)
end

function induce_acc(z::ComplexF64, ż::ComplexF64,
                           b::Vortex.Blob, blob_vel::ComplexF64)
    induce_acc(z, ż,
               Vortex.Point(Elements.position(b), circulation(b)),
               blob_vel)
end

function surface_pressure_inst(p::Plate, ṗ, ambient_sys, z_new, t, Δt, lesp, tesp)
    @get p (L, C, ss, α, dchebt!)
    @get ṗ (ċ, α̇)

    # Get Ċ from movement of existing vortex blobs (without vortex shedding)
    enforce_no_flow_through!(p, ṗ, ambient_sys, t)

    srcvel = self_induce_velocity(ambient_sys, t)
    induce_velocity!(srcvel, ambient_sys, p, t)

    targvel = fill(ċ, length(p))
    Ċ = zero(targvel)

    induce_acc!(Ċ, p.zs, targvel, ambient_sys, srcvel)

    z₊, z₋ = z_new
    point₊ = Vortex.Point(z₊, 1.0)
    point₋ = Vortex.Point(z₋, 1.0)

    Γ₊, Γ₋, ∂C₊, ∂C₋ = vorticity_flux!(p, point₊, point₋, t, lesp, tesp)

    n̂ = exp(-im*α)
    rmul!(Ċ, n̂)
    dchebt! * Ċ

    @. Ċ += (∂C₊ + ∂C₋)/Δt - im*α̇*C

    Ȧ = MappedVector(imag, Ċ, 1)

    # Plate is moving at a fixed velocity and angle of attack
    Ḃ₀ = Ḃ₁ = 0.0

    Γ̇ = [Γ₋/Δt + _bound_circulation(Ȧ, Ḃ₀, Ḃ₁, L, -(Γ₊ + Γ₋)/Δt, s) for s in p.ss]

    strength(p) .* (Chebyshev.firstkind(real.(C), ss) .- tangent(ċ, α)) .+ Γ̇
end
