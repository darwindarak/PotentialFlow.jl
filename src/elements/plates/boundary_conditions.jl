"""
    enforce_no_flow_through!(p::Plate, elements)

Update the plate, `p`, to enforce the no-flow-through condition given a collection of ambient vortex elements, `elements`.
"""
function enforce_no_flow_through!(p::Plate, elements)
    @get p (C, α)
    
    fill!(C, zero(Complex128))
    Vortex.induce_velocity!(C, p, elements)

    n̂ = exp(-im*α)
    scale!(C, n̂)

    chebyshev_transform!(C)

    p.Γ = -Vortex.circulation(elements)
    nothing
end


"""
    influence_on_plate!(∂A::Vector{Complex128}, plate::Plate, v)

Compute in-place the change in plate's Chebyshev coefficients `∂A` by a vortex element `v`
"""
function influence_on_plate!(∂A::Vector{Complex128}, plate::Plate, v)
    fill!(∂A, zero(Complex128))
    Vortex.induce_velocity!(∂A, plate, v)
    scale!(∂A, exp(-im*plate.α))
    chebyshev_transform!(∂A)
    nothing
end

function influence_on_plate(plate::Plate, v)
    ∂A = Vector{Complex128}(plate.N)
    influence_on_plate!(∂A, plate, v)
    return ∂A
end

"""
    vorticity_flux(p::Plate, v₁, v₂, 
                   lesp = 0.0, tesp = 0.0,
                   ∂A₁ = Vector{Complex128}(plate.N),
                   ∂A₂ = Vector{Complex128}(plate.N))

Return strengths of new vortex elements that satisfies edge suction parameters.
For a given edge, if the current suction parameter is less than the criticial suction parameter, then no vorticity is released.  If it is higher, however, vorticity will be released so that the suction parameter equals the critical value.

# Arguments:

- `p`: the plate
- `v₁, v₂`: the vortex elements (with unit circulation) that the vorticity flux is going into
- `lesp`, `tesp`: the critical leading and trailing edge suction parameters we want to enforce.  By default, both parameters are set to 0.0 to enforce the Kutta condition on both edges.  We can disable vortex shedding from an edge by setting the its critical suction parameter to `Inf`

# Returns:

- `Γ₁, Γ₂`: the strengths that the vortex element should have in order to satisfy the edge suction parameters
- `∂A₁, ∂A₂`: Chebyshev coefficients of the normal velocity induced by the vortex elements 
  Instead of running `enforce_bc!` with the new vortex elements, we can use this matrix to directly update the Chebyshev coefficients associated with the bound vortex sheet without recomputing all the velocities.
"""
function vorticity_flux(plate::Plate, v₁, v₂, lesp = 0.0, tesp = 0.0,
                        ∂A₁ = Vector{Complex128}(plate.N),
                        ∂A₂ = Vector{Complex128}(plate.N))

    @get plate (N, α, ċ, α̇, L, A, Γ)

    B₀ = imag(exp(-im*α)*ċ)
    B₁ = 0.5α̇*L

    b₊ = +2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L) 
    b₋ = -2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L) 

    influence_on_plate!(∂A₁, plate, v₁)
    influence_on_plate!(∂A₂, plate, v₂)

    Γ₁ = Vortex.circulation(v₁)
    Γ₂ = Vortex.circulation(v₂)

    A₁₊ =  2imag(∂A₁[1]) + imag(∂A₁[2]) - 2Γ₁/(π*L) 
    A₂₊ =  2imag(∂A₂[1]) + imag(∂A₂[2]) - 2Γ₂/(π*L) 
    A₁₋ = -2imag(∂A₁[1]) + imag(∂A₁[2]) - 2Γ₁/(π*L) 
    A₂₋ = -2imag(∂A₂[1]) + imag(∂A₂[2]) - 2Γ₂/(π*L) 

    if (abs2(lesp) > abs2(b₊)) && (abs2(tesp) ≤ abs2(b₋))
        K₁, K₂ = 0.0, (tesp - b₋)/A₂₋
    elseif (abs2(lesp) ≤ abs2(b₊)) && (abs2(tesp) > abs2(b₋))
        K₁, K₂ = (lesp - b₊)/A₁₊, 0.0
    elseif (abs2(lesp) > abs2(b₊)) && (abs2(tesp) > abs2(b₋))
        K₁ = K₂ = 0.0
    else
        b₊ = lesp - b₊
        b₋ = tesp - b₋

        detA = A₁₊*A₂₋ - A₂₊*A₁₋

        if detA == 0
            @show b₊, b₋
            @show A₁₊, A₁₋, A₂₊, A₂₋
            error("Cannot enforce suction parameters")
        end

        K₁ = (A₂₋*b₊ - A₂₊*b₋)/detA
        K₂ = (A₁₊*b₋ - A₁₋*b₊)/detA
    end

    return K₁*Γ₁, K₂*Γ₂, scale!(∂A₁, K₁), scale!(∂A₂, K₂)
end

function suction_parameters(plate)
    @get plate (N, α, ċ, α̇, L, A, Γ)

    B₀ = normal(ċ, α)
    B₁ = 0.5α̇*L

    b₊ = +2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L) 
    b₋ = -2(A[0] - B₀) + (A[1] - B₁) + 2Γ/(π*L) 
    return b₊, b₋
end                                  
