import ..Vortex: @property

@property begin
    signature = rate_of_impulse(plate, motion, elements::Source, velocities::Source)
    reduce = (+)
    stype = Complex128
end

function rate_of_impulse(plate, motion, element, velocity, ::Type{Vortex.Singleton})
    @get plate (c, α, L)
    @get motion (ċ, α̇)

    Γ = Vortex.circulation(element)
    z = 2exp(-im*α)*(Vortex.position(element) - c)/L
    ż = 2exp(-im*α)*(velocity - ċ)/L - im*α̇*z
    p̂ = imag(z) - im*real(√(z - 1)*√(z + 1))

    dp̂dt = -im*(ż + real((z/(√(z - 1)*√(z + 1)) - 1)*ż))

    0.5L*Γ*exp(im*α)*(dp̂dt + im*α̇*p̂)
end

"""
    force(plate, motion, elements, velocities, newelements = ())

Compute the force on `plate`, given its motion and the state of the ambient vorticity.

# Arguments

- `plate`: the plate
- `motion`: a structure that contains the velocity, acceleration, and angular velocity
  of the plate.
- `elements`: vortex elements representing the ambient vorticity
- `velocities`: the velocities of the vortex elements
- `newelements`: an optional argument listing vortex elements that are just added to the
  flow field (it can be an element that is contained in `elements`)
- `Δt`: this is only required if `newelements` is not empty, we assume that the new
  vortex elements are created over the span of `Δt`

# Returns

- `F`: the force excerted on the plate in complex coordinates
"""
function force(plate, motion, elements, velocities, newelements = (), Δt = 0.0)
    @get plate (c, α, L)
    @get motion (ċ, α̇, c̈)

    Ṗa = rate_of_impulse(plate, motion, elements, velocities)

    !isempty(newelements) && Δt == 0 && error("Δt should not be zero")

    for v in newelements
        Γ̇ = Vortex.circulation(v)/Δt
        z = 2exp(-im*α)*(Vortex.position(v) - c)/L
        p̂ = imag(z) - im*real(√(z - 1)*√(z + 1))

        Ṗa += exp(im*α)*0.5L*Γ̇*p̂
    end

    Ṗb = -im*π*exp(im*α)*(L/2)^2*( α̇*conj(exp(-im*α)*ċ)
                         - imag(exp(-im*α)*c̈)
                         )
    return -(Ṗa + Ṗb)
end
