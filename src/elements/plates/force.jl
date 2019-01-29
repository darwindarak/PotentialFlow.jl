import ..Properties: @property

"""
    rate_of_impulse(plate, motion, elements::Source, velocities::Source)

Compute the rate of change of impulse of a vortex element and its
image relative to a plate.

Note that this is not just the rate of impulse of the vortex element
itself, but also includes the rate of impulse of the bound vortex
sheet generated in response to the vortex element.
"""
@property begin
    signature = rate_of_impulse(plate, motion, elements::Source, velocities::Source)
    reduce = (+)
    stype = ComplexF64
end

function rate_of_impulse(plate, motion, element, velocity, ::Type{Singleton})
    @get plate (c, α, L)
    @get motion (ċ, α̇)

    Γ = circulation(element)
    z = 2exp(-im*α)*(position(element) - c)/L
    ż = 2exp(-im*α)*(velocity - ċ)/L - im*α̇*z
    p̂ = imag(z) - im*real(√(z - 1)*√(z + 1))

    dp̂dt = -im*(ż + real((z/(√(z - 1)*√(z + 1)) - 1)*ż))

    0.5L*Γ*exp(im*α)*(dp̂dt + im*α̇*p̂)
end

"""
    vortex_impulse(plate, elements::Source)

Compute the impulse of a vortex element and its image relative to a plate.

Note that this is not just the impulse of the vortex element
itself, but also includes the impulse of the bound vortex
sheet generated in response to the vortex element.
"""
@property begin
    signature = vortex_impulse(plate, elements::Source)
    reduce = (+)
    stype = ComplexF64
end

function vortex_impulse(plate, element, ::Type{Singleton})
    @get plate (c, α, L)

    Γ = circulation(element)
    z = 2exp(-im*α)*(position(element) - c)/L
    p̂ = imag(z) - im*real(√(z - 1)*√(z + 1))

    0.5L*Γ*exp(im*α)*p̂
end

"""
    rate_of_angular_impulse(plate, motion, elements::Source, velocities::Source)

Compute the rate of change of angular impulse of a vortex element (about the centroid)
and its image relative to a plate.

Note that this is not just the rate of angular impulse of the vortex element
itself, but also includes the rate of angular impulse of the bound vortex
sheet generated in response to the vortex element.
"""
@property begin
    signature = rate_of_angular_impulse(plate, motion, elements::Source, velocities::Source)
    reduce = (+)
    stype = ComplexF64
end

function rate_of_angular_impulse(plate, motion, element, velocity, ::Type{Singleton})
    @get plate (c, α, L)
    @get motion (ċ, α̇)

    Γ = circulation(element)
    z = 2exp(-im*α)*(position(element) - c)/L
    ż = 2exp(-im*α)*(velocity - ċ)/L - im*α̇*z
    #Π = -0.5abs(z)^2 -0.5real(z*(√(z - 1)*√(z + 1)-z))

    dΠdt = -real(conj(z)*ż + ((z^2-1/2)/(√(z - 1)*√(z + 1)) - z)*ż)

    0.25L^2*Γ*dΠdt
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

- `F`: the force exerted on the plate in complex coordinates
"""
function force(plate, motion, elements, velocities, newelements = (), Δt = 0.0)
    @get plate (c, α, L)
    @get motion (ċ, α̇, c̈)

    Ṗa = rate_of_impulse(plate, motion, elements, velocities)

    !isempty(newelements) && Δt == 0 && error("Δt should not be zero")

    for v in newelements
        Γ̇ = circulation(v)/Δt
        z = 2exp(-im*α)*(position(v) - c)/L
        p̂ = imag(z) - im*real(√(z - 1)*√(z + 1))

        Ṗa += exp(im*α)*0.5L*Γ̇*p̂
    end

    Ṗb = -im*π*exp(im*α)*(L/2)^2*( α̇*conj(exp(-im*α)*ċ)
                         - imag(exp(-im*α)*c̈)
                         )
    return -(Ṗa + Ṗb)
end

"""
    moment(plate, motion, elements, velocities, newelements = ())

Compute the moment on `plate` about its centroid, given its motion and the state of the ambient vorticity.

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

- `M`: the moment exerted on the plate (about the centroid)
"""
function moment(plate, motion, elements, velocities, newelements = (), Δt = 0.0)
    @get plate (c, α, L)
    @get motion (ċ, α̇, c̈, α̈)

    c̃̇ = ċ*exp(-im*α)

    Π̇a = rate_of_angular_impulse(plate, motion, elements, velocities)
    Π̇a += imag(conj(ċ)*vortex_impulse(plate,elements))

    !isempty(newelements) && Δt == 0 && error("Δt should not be zero")

    for v in newelements
        Γ̇ = circulation(v)/Δt
        z = 2exp(-im*α)*(position(v) - c)/L
        Π = -0.5abs(z)^2 -0.5real(z*(√(z - 1)*√(z + 1)-z))

        Π̇a += 0.25L^2*Γ̇*Π
    end


    Π̇b = 0.125π*(L/2)^4*α̈ + π*(L/2)^2*real(c̃̇)*imag(c̃̇)
    return -real.(Π̇a + Π̇b)
end
