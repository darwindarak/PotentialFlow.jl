function Vortex.induce_acc(z::Complex128, ż::Complex128,
                           p::Vortex.Point, point_vel::Complex128)
    -p.Γ*Vortex.Points.cauchy_kernel((z - p.z)^2)*conj(ż - point_vel)
end

function Vortex.induce_acc(z::Complex128, ż::Complex128,
                           b::Vortex.Blob, blob_vel::Complex128)
    Vortex.induce_acc(z, ż, Vortex.Point(b), blob_vel)
end

# We can use `_bound_circulation(A, B₀, B₁, L, Γ, s)` to compute
# the rate of change of Γ by supplying Ȧ, Ḃ₀, Ḃ₁, Γ̇
#
# Ȧ can be computed

# fill!(C, zero(Complex128))
# Vortex.induce_velocity!(C, p, elements)
#
# n̂ = exp(-im*α)
# scale!(C, n̂)
#
# dchebt! * C
#
# p.Γ = -Vortex.circulation(elements)
# p.B₀ = normal(ċ, α)
# p.B₁ = 0.5α̇*L

# d/dt (exp(-im*α)*U) = -im*α̇*exp(-im*α)*U + exp(-im*α)*U̇
#
# -im*α̇*exp(-im*α)*U = -im*α̇*(Σ Cₙ Tₙ)
# exp(-im*α)*U̇ = Σ K̇ₙ Tₙ
#
# d/dt (exp(-im*α)*U) = Σ (K̇ₙ - im*α̇*Cₙ) Tₙ = Σ Ċₙ Tₙ

function surface_pressure_inst(p::Plate, ṗ, blobs, Δt, lesp, tesp)
    @get p (L, C, ss, α, dchebt!)
    @get ṗ (ċ, α̇)

    # Get Ċ from movement of existing vortex blobs (without vortex shedding)
    enforce_no_flow_through!(p, ṗ, blobs)

    srcvel = Vortex.self_induce_velocity(blobs)
    Vortex.induce_velocity!(srcvel, blobs, p)

    targvel = fill(ċ, length(p))
    Ċ = zeros(targvel)

    Vortex.induce_acc!(Ċ, p.zs, targvel, blobs, srcvel)

    δ = blobs[1].δ
    z₊ = (blobs[end-1].z + 2p.zs[end])/3
    z₋ = (blobs[end].z + 2p.zs[1])/3
    blob₊ = Vortex.Blob(z₊, 1.0, δ)
    blob₋ = Vortex.Blob(z₋, 1.0, δ)

    Γ₊, Γ₋, ∂C₊, ∂C₋ = vorticity_flux!(p, blob₊, blob₋, lesp, tesp)

    # Γ̇₊, Γ̇₋, ∂Ċ₊, ∂Ċ₋ = circulation_shed_rate(p, ṗ, blobs, Δt, lesp, tesp)

    n̂ = exp(-im*α)
    scale!(Ċ, n̂)
    dchebt! * Ċ

    @. Ċ += (∂C₊ + ∂C₋)/Δt - im*α̇*C

    Ȧ = MappedVector(imag, Ċ, 1)

    # Plate is moving at a fixed velocity and angle of attack
    Ḃ₀ = Ḃ₁ = 0.0

    Γ̇ = [Γ₋/Δt + _bound_circulation(Ȧ, Ḃ₀, Ḃ₁, L, -(Γ₊ + Γ₋)/Δt, s) for s in p.ss]

    strength(p) .* (Chebyshev.firstkind(real.(C), ss) .- tangent(ċ, α)) .+ Γ̇
end

function circulation_shed_rate(plate, motion, blobs, Δt, lesp = 0.0, tesp = 0.0)
    δ = blobs[1].δ
    z₊ = (blobs[end-1].z + 2plate.zs[end])/3
    z₋ = (blobs[end].z + 2plate.zs[1])/3

    blob₊ = Vortex.Blob(z₊, 1.0, δ)
    blob₋ = Vortex.Blob(z₋, 1.0, δ)
    Vortex.Plates.enforce_no_flow_through!(plate, motion, blobs)

    Vortex.Plates.vorticity_flux!(plate, blob₊, blob₋, lesp, tesp) ./ Δt
end
