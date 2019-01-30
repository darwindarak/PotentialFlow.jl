# PotentialFlow

*a scaffolding for building 2D inviscid models*

| Documentation | Build Status |
|:---:|:---:|
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://darwindarak.github.com/PotentialFlow.jl) [![Binder](https://beta.mybinder.org/badge.svg)](https://beta.mybinder.org/v2/gh/darwindarak/PotentialFlow.jl/binder?filepath=binder%2Fnotebooks) | [![Build Status](https://img.shields.io/travis/darwindarak/PotentialFlow.jl/master.svg?label=linux)](https://travis-ci.org/darwindarak/PotentialFlow.jl) [![Build status](https://img.shields.io/appveyor/ci/darwindarak/PotentialFlow-jl/master.svg?label=windows)](https://ci.appveyor.com/project/darwindarak/potentialflow-jl/branch/master) [![codecov](https://codecov.io/gh/darwindarak/PotentialFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/darwindarak/PotentialFlow.jl) |

## Installation

PotentialFlow can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia-repl
pkg> add PotentialFlow
```

## Basic Usage

Let's create a point vortex and a point source and probe their associated velocity field:
```julia-repl
julia> using PotentialFlow

julia> t = 0.0
0.0

julia> vortex = Vortex.Point(1.0 + 1.0im, π)
Vortex.Point(1.0 + 1.0im, 3.141592653589793)

julia> source = Source.Point(1.0 - 1.0im, -π)
Source.Point(1.0 - 1.0im, 3.141592653589793)

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
```julia
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

x = range(-2, 1, length=100)
y = range(-0.5, 0.5, length=100)

streamlines(x, y, (plate, freestream), legend = false, colorbar = false)
plot!(plate, linewidth = 2, ratio = 1, size = (600, 300))
```
![Flat plate in freestream](https://darwindarak.github.io/PotentialFlow.jl/latest/translating_plate.svg)

More examples can be found in the [documentation](https://darwindarak.github.io/PotentialFlow.jl) and the [Jupyter notebooks](https://github.com/darwindarak/PotentialFlow.jl/tree/binder/binder/notebooks).
You can also run the notebooks directly in your browser [here](https://beta.mybinder.org/v2/gh/darwindarak/PotentialFlow.jl/binder?filepath=binder%2Fnotebooks).
