# Plate Motions

```@meta
CurrentModule = Plates.RigidBodyMotions
DocTestSetup  = quote
    using PotentialFlow
    using Random
    Random.seed!(1)
end
```

The motion of a plate is specified through two data types:
- [`RigidBodyMotion`](@ref) is the type that should be used to represent the
  plate's velocity.  For example, in `advect!(plate₊, plate₋,
  platevel, Δt)`, `platevel` is of type `RigidBodyMotion.` It contains the most
  current values `(ċ, c̈, α̇)` (the plate's centroid velocity and
  acceleration, and angular velocity, respectively), as well as a
  [`Kinematics`](@ref) type.

- [`Kinematics`](@ref) is an abstract type representing a function
  that takes in a time and returns `(ċ, c̈, α̇)`

## Motion

By default, `RigidBodyMotion` assumes a constant translational and angular velocity.
For example,
```jldoctest constant
julia> motion = Plates.RigidBodyMotion(1.0im, π/2)
Rigid Body Motion:
  ċ = 0.0 + 1.0im
  c̈ = 0.0 + 0.0im
  α̇ = 1.57
  α̈ = 0.0
  Constant (ċ = 0.0 + 1.0im, α̇ = 1.5707963267948966)
```
Here, `Constant` is a subtype of [`Kinematics`](@ref) that returns the same `(ċ, c̈, α̇)` triple at all times
```jldoctest constant
julia> motion.kin.([0.0, 1.0, 2.0])
3-element Array{Tuple{Complex{Float64},Complex{Float64},Float64,Complex{Float64}},1}:
 (0.0 + 1.0im, 0.0 + 0.0im, 1.5707963267948966, 0.0 + 0.0im)
 (0.0 + 1.0im, 0.0 + 0.0im, 1.5707963267948966, 0.0 + 0.0im)
 (0.0 + 1.0im, 0.0 + 0.0im, 1.5707963267948966, 0.0 + 0.0im)

```
Calling `Plates.RigidBodyMotion(1.0im, π/2)` is equivalent doing
```jldoctest
kin = Plates.RigidBodyMotions.Constant(1.0im, π/2)
motion = Plates.RigidBodyMotion(1.0im, 0.0im, π/2, 0.0, kin)

# output

Rigid Body Motion:
  ċ = 0.0 + 1.0im
  c̈ = 0.0 + 0.0im
  α̇ = 1.57
  α̈ = 0.0
  Constant (ċ = 0.0 + 1.0im, α̇ = 1.5707963267948966)
```
The next section describes how to construct more interesting kinematics.

## Kinematics

The [`Kinematics`](@ref) type is just an abstract type for functions
that take in time and return the `(ċ, c̈, α̇)` triple.  Let's create a
`MyMotion` type that describes a horizontally translating plate that
also sinusoidally pitches about its centroid.
```jldoctest sinusoidal
import PotentialFlow.Plates.RigidBodyMotions: Kinematics

struct MyMotion <: Kinematics
    U₀::ComplexF64
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
 (1.0 + 0.0im, 0.0 + 0.0im, 0.0)
 (1.0 + 0.0im, 0.0 + 0.0im, 0.7071067811865475)
 (1.0 + 0.0im, 0.0 + 0.0im, 1.0)

```

## Profiles

To make defining complex kinematics a little eaiser, the library also
provides a [`Profile`](@ref) type, an abstract type for
real-valued functions of time.
Before going into how to define new profiles, we'll first show an
example of why we might want to represent functions as a type.
We start off with a predefined profile, a smooth ramp:
```@example ramp
using Plots
using PotentialFlow.Plates.RigidBodyMotions

ramp = RigidBodyMotions.EldredgeRamp(6)

T = range(-1, 4, length=200)
plot(T, ramp.(T), xlabel = "t", ylabel="Smoothed Ramp",
     legend = :none, linewidth = 2)

savefig("ramp.svg"); nothing # hide
```
```@raw html
<object data="manual/ramp.svg" type="image/svg+xml"></object>
```
Now suppose we want to scale the ramp and shift it
```@example ramp
shifted_ramp = -(ramp >> 2)

plot(T, shifted_ramp.(T), xlabel = "t", ylabel="Smoothed Ramp",
     legend = :none, linewidth = 2, size=(600,300))
savefig("shifted_ramp.svg"); nothing # hide
```
```@raw html
<object data="manual/shifted_ramp.svg" type="image/svg+xml"></object>
```
then take its derivative
```@example ramp
ddt_ramp = d_dt(shifted_ramp)

plot(T, ddt_ramp.(T), xlabel = "t", ylabel="Smoothed Ramp",
     legend = :none, linewidth = 2, size = (600, 200))
savefig("ddt_ramp.svg"); nothing # hide
```
```@raw html
<object data="manual/ddt_ramp.svg" type="image/svg+xml"></object>
```
We see that wrapping these functions in a type allows us to operate on
them as if they values, making it easier to compose multiple motions together:
```@example ramp
ps_ramp = RigidBodyMotions.ColoniusRamp(5)
composed_ramp = ramp - (ps_ramp >> 2)

plot(T, composed_ramp.(T), xlabel = "t", ylabel="Smoothed Ramp",
     legend = :none, linewidth = 2, size = (600, 300))
savefig("composed_ramp.svg"); nothing # hide
```
```@raw html
<object data="manual/composed_ramp.svg" type="image/svg+xml"></object>
```

### Defining a profile

Defining a profile is done in two steps:

1. Create a subtype of `RigidBodyMotions.Profile` that contains the relavant parameters, e.g.
2. Add a method on the type (see [Function like objects](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1))

For example,
```@example custom_profile
using PotentialFlow.Plates.RigidBodyMotions

struct Sinusoid <: RigidBodyMotions.Profile
    ω::Float64
end

(s::Sinusoid)(t) = sin(s.ω*t)
```
which can then be used as follows:
```@example custom_profile

T = range(-6, 6, length = 200)

s = Sinusoid(2.0)
c = d_dt(2s >> 0.5)

using Plots
plot(T, [s.(T) c.(T)], xlabel = "t", color = ["#00BFFF" "#D4CA3A"],
     legend = :none, linewidth = 2)
savefig("custom_profile.svg"); nothing # hide
```
```@raw html
<object data="manual/custom_profile.svg" type="image/svg+xml"></object>
```

## Function Documentation

```@autodocs
Modules = [RigidBodyMotions]
Order   = [:type, :function]
```

## Index

```@index
Pages   = ["motions.md"]
```
