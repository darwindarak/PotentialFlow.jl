# VortexModel

*a scaffolding for building vortex models*

The main goal of this library is to remove as much boilerplate code as possible from vortex modeling codes.
The core operation in vortex models is simulating the dynamics of various interacting vortex elements.
In general, the simulation comes down to computing the velocities of the vortex elements then applying some time-marching scheme to evolve the system forward in time.
With this in mind, we want to construct a library that makes it

- easy to define new vortex types and behaviors
- straightforward for users to set up a system the vortex elements
- intuitive to probe the state of any vortex element in the system
- easy to define new time-marching schemes to fit the users needs

## Installation

This package requires Julia `0.6-` and above.
It is not a registered package, so it should be installed with:
```julia
julia> Pkg.clone("git@github.com:darwindarak/VortexModel.jl.git")
```
Since it is still under heavy development, you should run
```julia
julia> Pkg.test("VortexModel")
```
to make sure things are working as intended.

The plots in this documentation are generated using [PyPlot.jl](github.com/JuliaPy/PyPlot.jl).
You might want to install that too to follow the examples in the [getting started guide](@ref getting-started) or the [Jupyter notebooks](https://github.com/darwindarak/VortexModel.jl/tree/master/examples).
