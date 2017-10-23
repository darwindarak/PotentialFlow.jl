# PotentialFlow

*a scaffolding for building 2D inviscid models*

The main goal of this library is to remove as much boilerplate code as possible from inviscid modeling codes.
The core operation in these models is simulating the dynamics of various interacting potential flow elements.
In general, the simulation comes down to computing the velocities of the elements then applying some time-marching scheme to evolve the system forward in time.
With this in mind, we want to construct a library that makes it

- easy to define new flow elements and behaviors
- straightforward for users to set up a system of elements
- intuitive to probe the state of any element in the system
- easy to define new time-marching schemes to fit the users needs

## Installation

This package requires Julia `0.6-` and above.
To install, simply run
```julia
julia> Pkg.add("PotentialFlow")
```
in the Julia REPL.
Since this package is still under heavy development, you should run
```julia
julia> Pkg.test("PotentialFlow") # might take some time
```
to make sure things are working as intended and
```julia
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that too to follow the examples in the [getting started guide](@ref getting-started) or the [Jupyter notebooks](https://github.com/darwindarak/PotentialFlow.jl/tree/master/binder/notebooks).
The notebooks can also be run online [here](https://beta.mybinder.org/v2/gh/darwindarak/PotentialFlow.jl/master?filepath=binder%2Fnotebooks).

## Basic Usage

Let's create a point vortex and a point source and probe their associated velocity field:
```jldoctest
julia> using PotentialFlow

julia> t = 0.0
0.0

julia> vortex = Vortex.Point(1.0 + 1.0im, π)
Vortex.Point(1.0 + 1.0im, 3.141592653589793)

julia> source = Source.Point(1.0 - 1.0im, -π)
Source.Point(1.0 - 1.0im, -3.141592653589793)

julia> induce_velocity(0.0im, vortex, t)
0.25 - 0.25im

julia> induce_velocity(source, vortex, t)
0.25 - 0.0im

julia> induce_velocity(0.0im, (vortex, source), t)
0.5 - 0.5im

julia> induce_velocity([0.0im, 1.0im, 1.0], (vortex, source), t)
3-element Array{Complex{Float64},1}:
 0.5-0.5im
 0.1-0.7im
 0.5-0.5im
```
Note the all positions and velocities are given in complex coordiantes.

Now let's move on to something more interesting.
We'll create a stationary flat plate (bound vortex sheet) and place it in a freestream.
In order to enforce the Kutta condition, we also place a starting vortex at `-Inf`.
```@example plate-kutta
using PotentialFlow
using Plots

c₀ = 0.0im # initial centroid position
α = π/9    # angle of attack
L = 1.0    # chord length
N = 128    # number of discretization points

ċ = 0.0    # translation velocity
α̇ = 0.0    # rate of rotation
t = 0.0    # current time

freestream = Freestream(-1.0)

plate = Plate(N, L, c₀, α)
motion = Plates.RigidBodyMotion(ċ, α̇)
Plates.enforce_no_flow_through!(plate, motion, freestream, 0.0)

# We now want to determine the strength of the starting vortex
# to satisfy the Kutta condition at the trailing edge of the plate
_, Γ = Plates.vorticity_flux!(plate, (), Vortex.Point(-Inf, 1.0), t, Inf, 0);
starting_vortex = Vortex.Point(-Inf, Γ)

# Plot some streamlines

x = linspace(-2, 1, 100)
y = linspace(-0.5, 0.5, 100)

streamlines(x, y, (plate, freestream), legend = false, colorbar = false)
plot!(plate, linewidth = 2, ratio = 1, size = (600, 300))
savefig("translating_plate.svg") # hide
nothing # hide
```
![Flat plate in freestream](translating_plate.svg)
