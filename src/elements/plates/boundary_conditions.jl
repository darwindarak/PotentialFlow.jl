"""
    enforce_no_flow_through!(p::Plate, motion, elements, t)

Update the plate, `p`, to enforce the no-flow-through condition given ambient vortex elements, `elements`, and while moving with kinematics specified by `motion`.

# Example

```jldoctest
julia> plate = Plate(128, 2.0, 0.0, π/3)
Plate: N = 128, L = 2.0, c = 0.0 + 0.0im, α = 60.0ᵒ
       LESP = 0.0, TESP = 0.0

julia> motion = Plates.RigidBodyMotion(1.0, 0.0);

julia> point = Vortex.Point(0.0 + 2im, 1.0);

julia> Plates.enforce_no_flow_through!(plate, motion, point, 0.0)

julia> plate
Plate: N = 128, L = 2.0, c = 0.0 + 0.0im, α = 60.0ᵒ
       LESP = 1.27, TESP = -1.93
```
"""
function enforce_no_flow_through!(p::Plate, ṗ, elements, t)
    @get p (L, C, α, dchebt!)
    @get ṗ (ċ, α̇)

    fill!(C, zero(ComplexF64))
    induce_velocity!(C, p, elements, t)

    n̂ = exp(-im*α)
    rmul!(C, n̂)

    dchebt! * C

    p.Γ = -circulation(elements)
    p.B₀ = normal(ċ, α)
    p.B₁ = 0.5α̇*L

    nothing
end


"""
    influence_on_plate!(∂A::Vector{ComplexF64}, plate::Plate, v, t)

Compute in-place the change in plate's Chebyshev coefficients `∂A` by a vortex element `v`
"""
function influence_on_plate!(∂A::Vector{ComplexF64}, plate::Plate, v, t)
    fill!(∂A, zero(ComplexF64))
    induce_velocity!(∂A, plate, v, t)
    rmul!(∂A, exp(-im*plate.α))
    plate.dchebt! * ∂A
    nothing
end

function influence_on_plate(plate::Plate, v, t)
    ∂A = Vector{ComplexF64}(undef, plate.N)
    influence_on_plate!(∂A, plate, v, t)
    return ∂A
end

"""
    vorticity_flux(p::Plate, v₁, v₂,
                   lesp = 0.0, tesp = 0.0,
                   ∂C₁ = Vector{ComplexF64}(undef, plate.N),
                   ∂C₂ = Vector{ComplexF64}(undef, plate.N))

Return strengths of new vortex elements that satisfies edge suction parameters.
For a given edge, if the current suction parameter is less than the criticial suction parameter, then no vorticity is released.  If it is higher, however, vorticity will be released so that the suction parameter equals the critical value.

# Arguments

- `p`: the plate
- `v₁, v₂`: the vortex elements (with unit circulation) that the vorticity flux is going into
- `lesp`, `tesp`: the critical leading and trailing edge suction parameters we want to enforce.  By default, both parameters are set to 0.0 to enforce the Kutta condition on both edges.  We can disable vortex shedding from an edge by setting the its critical suction parameter to `Inf`

# Returns

- `Γ₁, Γ₂`: the strengths that the vortex element should have in order to satisfy the edge suction parameters
- `∂C₁, ∂C₂`: Chebyshev coefficients of the normal velocity induced by the vortex elements
  Instead of running `enforce_bc!` with the new vortex elements, we can use this matrix to directly update the Chebyshev coefficients associated with the bound vortex sheet without recomputing all the velocities.

# Example

Enforcing the trailing edge Kutta condition with an point vortex at negative infinity:

```jldoctest
julia> plate = Plate(128, 2.0, 0.0, π/6)
Plate: N = 128, L = 2.0, c = 0.0 + 0.0im, α = 30.0ᵒ
       LESP = 0.0, TESP = 0.0

julia> motion = Plates.RigidBodyMotion(1.0, 0.0);

julia> Plates.enforce_no_flow_through!(plate, motion, (), 0.0)

julia> point = Vortex.Point(-Inf, 1.0);

julia> _, Γ, _, _ = Plates.vorticity_flux(plate, (), point, 0.0, Inf);

julia> Γ # should equal -πULsin(α) = -π
-3.1415926535897927
```
"""
function vorticity_flux(plate::Plate, v₁, v₂, t, lesp = 0.0, tesp = 0.0,
                        ∂C₁ = Vector{ComplexF64}(undef,plate.N),
                        ∂C₂ = Vector{ComplexF64}(undef,plate.N))

    @get plate (N, α, B₀, B₁, L, A, Γ)

    b₊ = +2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L)
    b₋ = -2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L)

    influence_on_plate!(∂C₁, plate, v₁, t)
    influence_on_plate!(∂C₂, plate, v₂, t)

    Γ₁ = circulation(v₁)
    Γ₂ = circulation(v₂)

    A₁₊ =  2imag(∂C₁[1]) + imag(∂C₁[2]) - 2Γ₁/(π*L)
    A₂₊ =  2imag(∂C₂[1]) + imag(∂C₂[2]) - 2Γ₂/(π*L)
    A₁₋ = -2imag(∂C₁[1]) + imag(∂C₁[2]) - 2Γ₁/(π*L)
    A₂₋ = -2imag(∂C₂[1]) + imag(∂C₂[2]) - 2Γ₂/(π*L)

    if (abs2(lesp) > abs2(b₊)) && (abs2(tesp) ≤ abs2(b₋))
        K₁, K₂ = 0.0, (sign(b₋)*tesp - b₋)/A₂₋
    elseif (abs2(lesp) ≤ abs2(b₊)) && (abs2(tesp) > abs2(b₋))
        K₁, K₂ = (sign(b₊)*lesp - b₊)/A₁₊, 0.0
    elseif (abs2(lesp) > abs2(b₊)) && (abs2(tesp) > abs2(b₋))
        K₁ = K₂ = 0.0
    else
        b₊ = sign(b₊)*lesp - b₊
        b₋ = sign(b₋)*tesp - b₋

        detA = A₁₊*A₂₋ - A₂₊*A₁₋

        @assert (detA != 0) "Cannot enforce suction parameters"

        K₁ = (A₂₋*b₊ - A₂₊*b₋)/detA
        K₂ = (A₁₊*b₋ - A₁₋*b₊)/detA
    end

    return K₁*Γ₁, K₂*Γ₂, rmul!(∂C₁, K₁), rmul!(∂C₂, K₂)
end

"""
    vorticity_flux!(p::Plate, v₁, v₂,
                    lesp = 0.0, tesp = 0.0,
                    ∂C₁ = Vector{ComplexF64}(undef,plate.N),
                    ∂C₂ = Vector{ComplexF64}(undef,plate.N))

In-place version of [`vorticity_flux`](@ref), except instead of just
returning the possible changes in plate Chebyshev coefficients, we
modify `plate.C` with those changes so that no-flow-through is
enforced in the presence of `v₁` and `v₂` with strengths that satisfy
the suction parameters.
"""
function vorticity_flux!(plate::Plate, v₁, v₂, t, lesp = 0.0, tesp = 0.0,
                        ∂C₁ = Vector{ComplexF64}(undef, plate.N),
                        ∂C₂ = Vector{ComplexF64}(undef, plate.N))
    Γ₁, Γ₂, _, _ = vorticity_flux(plate, v₁, v₂, t, lesp, tesp, ∂C₁, ∂C₂)
    @. plate.C += ∂C₁ + ∂C₂
    plate.Γ -= Γ₁ + Γ₂

    return Γ₁, Γ₂, ∂C₁, ∂C₂
end

function suction_parameters(plate)
    @get plate (N, α, B₀, B₁, L, A, Γ)

    b₊ = +2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L)
    b₋ = -2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L)
    return b₊, b₋
end
