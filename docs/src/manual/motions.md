# Plate Motions

```@meta
CurrentModule = Vortex.Plates.Motions
DocTestSetup  = quote
    import VortexModel.Vortex.Plates: Motions
    using .Motions
    srand(1)
end
```

The motion of a plate is specified through two data types:
- [`Motion`](@ref) is the type that should be used to represent the
  plate's velocity.  For example, in `advect!(plate₊, plate₋,
  platevel, Δt)`, `platevel` is of type `Motion.` It contains the most
  current values `(ċ, c̈, α̇)` (the plate's centroid velocity and
  acceleration, and angular velocity, respectively), as well as a
  [`Kinematics`](@ref) type.

- [`Kinematics`](@ref) is an abstract type representing a function
  that takes in a time and returns `(ċ, c̈, α̇)`

## Motion

By default, `Motion` assumes a constant translational and angular velocity.
For example,
```jldoctest constant
julia> motion = Motion(1.0im, π/2)
Plate Motion:
  ċ = 0.0 + 1.0im
  c̈ = 0.0 + 0.0im
  α̇ = 1.57
  Constant (ċ = 0.0 + 1.0im, α̇ = 1.5707963267948966)
```
Here, `Constant` is a subtype of [`Kinematics`](@ref) that returns the same `(ċ, c̈, α̇)` triple at all times
```jldoctest constant
julia> motion.kin.([0.0, 1.0, 2.0])
3-element Array{Tuple{Complex{Float64},Complex{Float64},Float64},1}:
 (0.0+1.0im, 0.0+0.0im, 1.5708)
 (0.0+1.0im, 0.0+0.0im, 1.5708)
 (0.0+1.0im, 0.0+0.0im, 1.5708)
```
Calling `Motion(1.0im, π/2)` is equivalent doing
```jldoctest
kin = Motions.Constant(1.0im, π/2)
motion = Motion(1.0im, 0.0im, π/2, kin)

# output

Plate Motion:
  ċ = 0.0 + 1.0im
  c̈ = 0.0 + 0.0im
  α̇ = 1.57
  Constant (ċ = 0.0 + 1.0im, α̇ = 1.5707963267948966)
```
The next section describes how to construct more interesting kinematics.

## Kinematics

The [`Kinematics`](@ref) type is just an abstract type for functions
that take in time and return the `(ċ, c̈, α̇)` triple.  Let's create a
`MyMotion` type that describes a horizontally translating plate that
also sinusoidally pitches about its centroid.
```jldoctest sinusoidal
struct MyMotion <: Kinematics
    U₀::Complex128
    ω::Float64
end

(m::MyMotion)(t) = (m.U₀, 0.0im, sin(m.ω*t))

sinusoid = MyMotion(1.0, π/4)

# output

MyMotion(1.0 + 0.0im, 0.7853981633974483)
```
We can then evaluate `sinusoid` at different times
```jldoctest sinusoidal
julia> sinusoid.([0.0, 1.0, 2.0])
3-element Array{Tuple{Complex{Float64},Complex{Float64},Float64},1}:
 (1.0+0.0im, 0.0+0.0im, 0.0)
 (1.0+0.0im, 0.0+0.0im, 0.707107)
 (1.0+0.0im, 0.0+0.0im, 1.0)
```

## Profiles

```@setup ramp
import VortexModel.Vortex.Plates: Motions
import .Motions: Kinematics, Motion, d_dt
using Gadfly
srand(1)
```
To make defining complex kinematics a little eaiser, the library also
provides a [`Motions.Profile`](@ref) type, an abstract type for
real-valued functions of time.
Before going into how to define new profiles, we'll first show an
example of why we might want to represent functions as a type.
We start off with a predefined profile, a smooth ramp:
```@example ramp
ramp = Motions.EldredgeRamp(6)

T = linspace(-1, 4, 200)
plot(x = T, y = ramp.(T), Geom.line, Guide.xlabel("t"))
draw(SVGJS("ramp.svg", 6inch, 4inch), ans); nothing # hide
```
```@raw html
<object data="ramp.svg" type="image/svg+xml"></object>
```
Now suppose we want to scale the ramp and shift it
```@example ramp
shifted_ramp = -(ramp >> 2)

plot(x = T, y = shifted_ramp.(T), Geom.line, Guide.xlabel("t"))
draw(SVGJS("shifted_ramp.svg", 6inch, 4inch), ans); nothing # hide
```
```@raw html
<object data="shifted_ramp.svg" type="image/svg+xml"></object>
```
then take its derivative
```@example ramp
ddt_ramp = d_dt(shifted_ramp)

plot(x = T, y = ddt_ramp.(T), Geom.line, Guide.xlabel("t"), Guide.ylabel("ẏ"))
draw(SVGJS("ddt_ramp.svg", 6inch, 4inch), ans); nothing # hide
```
```@raw html
<object data="ddt_ramp.svg" type="image/svg+xml"></object>
```
We see that wrapping these functions in a type allows us to operate on
them as if they values, making it easier to compose multiple motions together:
```@example ramp
ps_ramp = Motions.ColoniusRamp(5)
composed_ramp = ramp - (ps_ramp >> 2)

plot(x = T, y = composed_ramp.(T), Geom.line, Guide.xlabel("t"), Guide.ylabel("y"))
draw(SVGJS("composed_ramp.svg", 6inch, 4inch), ans); nothing # hide
```
```@raw html
<object data="composed_ramp.svg" type="image/svg+xml"></object>
```

### Defining a profile

```@setup custom_profile
import VortexModel.Vortex.Plates: Motions
import .Motions: Kinematics, Motion, d_dt
using Gadfly
srand(1)
```

Defining a profile is done in two steps:

1. Create a subtype of `Motions.Profile` that contains the relavant parameters, e.g.
2. Add a method on the type (see [Function like objects](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1))

For example,
```@example custom_profile
struct Sinusoid <: Motions.Profile
    ω::Float64
end

(s::Sinusoid)(t) = sin(s.ω*t)
```
which can then be used as follows:
```@example custom_profile
T = linspace(-6, 6, 200)

s = Sinusoid(2.0)
c = d_dt(2s >> 0.5)

plot(layer(x = T, y = s.(T), Geom.line, Theme(default_color = "#00BFFF")),
     layer(x = T, y = c.(T), Geom.line, Theme(default_color = "#D4CA3A")))
draw(SVGJS("custom_profile.svg", 6inch, 4inch), ans); nothing # hide
```
```@raw html
<object data="custom_profile.svg" type="image/svg+xml"></object>
```

## Function Documentation

```@autodocs
Modules = [Motions]
Order   = [:type, :function]
```

## Index

```@index
Pages   = ["motions.md"]
```
