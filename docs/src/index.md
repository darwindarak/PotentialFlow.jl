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
It is not a registered package, so it should be installed with:
```julia-repl
julia> Pkg.clone("git@github.com:darwindarak/PotentialFlow.jl.git")
```
Since it is still under heavy development, you should run
```julia-repl
julia> Pkg.test("PotentialFlow") # might take some time
```
to make sure things are working as intended and
```julia-repl
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that too to follow the examples in the [getting started guide](@ref getting-started) or the [Jupyter notebooks](https://github.com/darwindarak/PotentialFlow.jl/tree/master/examples).
