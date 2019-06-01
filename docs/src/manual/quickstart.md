# [Getting Started](@id getting-started)

This getting started guide will introduce the main components of **PotentialFlow.jl**.
The code examples here should be directly copy-paste-able into the Julia REPL (even with the `julia>` prompt and sample results).

```@meta
DocTestSetup = quote
    using Random
    Random.seed!(1)
end
```

## Creating Flow Elements

We start by importing the library and creating a single point vortex with unit circulation located at (1,1):
```jldoctest quickstart
julia> using PotentialFlow

julia> p = Vortex.Point( 1.0 + 1.0im, 1.0 )
Vortex.Point(1.0 + 1.0im, 1.0)
```
By convention, the arguments for element constructors are position(s), circulation/strength(s), followed by any type specific parameters.
For example, a vortex blob at the same location as `p` with a blob radius of 0.1 is created with
```jldoctest quickstart
julia> Vortex.Blob(1.0 + 1.0im, 1.0, 0.1)
Vortex.Blob(1.0 + 1.0im, 1.0, 0.1)
```

We can use Julia's [vectorized dot syntax](https://docs.julialang.org/en/latest/manual/functions.html#man-vectorized-1) to construct whole arrays of elements.
For example, here we create five point vortices and five point sources:
```jldoctest quickstart
julia> N = 5;

julia> zs = Complex.(randn(N), randn(N));

julia> vortices = Vortex.Point.(zs .+ 1.5, rand(N))
5-element Array{PotentialFlow.Points.Point{Float64},1}:
 Vortex.Point(1.7972879845354617 + 0.31111133849833383im, 0.42471785049513144)
 Vortex.Point(1.882395967790608 + 2.2950878238373105im, 0.773223048457377)
 Vortex.Point(0.9023655232717689 - 2.2670863488005306im, 0.2811902322857298)
 Vortex.Point(1.4895547553626243 + 0.5299655761667461im, 0.20947237319807077)
 Vortex.Point(0.660973145611236 + 0.43142152642291204im, 0.25137920979222494)

julia> sources = Source.Point.(zs .- 1.5, rand(N))
5-element Array{PotentialFlow.Points.Point{Complex{Float64}},1}:
 Source.Point(-1.2027120154645383 + 0.31111133849833383im, 0.02037486871266725)
 Source.Point(-1.117604032209392 + 2.2950878238373105im, 0.2877015122756894)
 Source.Point(-2.0976344767282313 - 2.2670863488005306im, 0.859512136087661)
 Source.Point(-1.5104452446373757 + 0.5299655761667461im, 0.07695088688120899)
 Source.Point(-2.339026854388764 + 0.43142152642291204im, 0.6403962459899388)

```

We can mix different vortex types together by grouping them in tuples.
For example, a collection of vortex elements consisting of the point vortices and vortex blobs created earlier can be grouped together with:
```jldoctest quickstart
julia> sys = (vortices, sources);
```

!!! note
    The Unicode characters used in the examples can be entered in the Julia REPL (and most text editors with the appropriate plugins) via [tab completion.](https://docs.julialang.org/en/latest/manual/unicode-input.html#Unicode-Input-1).  For example:
    - Γ: `\Gamma<TAB>`
    - Δ: `\Delta<TAB>`
    - ẋ: `x\dot<TAB>`
    - 🌀: `\:cyclone:<TAB>`


We can access properties of any vortex element by directly accessing its fields, for example:
```jldoctest quickstart
julia> p.z
1.0 + 1.0im

```
However, it is better practice to use accessor methods, such as:
```jldoctest quickstart
julia> Elements.position(p)
1.0 + 1.0im

```
since not all element types store their position in a `z` field but they are all required to implement a `Elements.position` method (also see `Elements.impulse` and `Elements.position`).
These accessor methods, combined with the dot syntax, also make it easier to work with properties of arrays and tuples of vortex elements.

```jldoctest quickstart
julia> Elements.circulation(vortices)
1.939982714228534

julia> Elements.circulation(sources)
0.0

julia> Elements.circulation(sys)
1.939982714228534

julia> Elements.circulation.(vortices)
5-element Array{Float64,1}:
 0.42471785049513144
 0.773223048457377
 0.2811902322857298
 0.20947237319807077
 0.25137920979222494

julia> Elements.position.(sources)
5-element Array{Complex{Float64},1}:
 -1.2027120154645383 + 0.31111133849833383im
  -1.117604032209392 + 2.2950878238373105im
 -2.0976344767282313 - 2.2670863488005306im
 -1.5104452446373757 + 0.5299655761667461im
  -2.339026854388764 + 0.43142152642291204im

```

## Computing Velocities

Now that we can create potential flow elements, we want to add in some dynamics.
The key functions for this are the `induce_velocity` and `induce_velocity!` pair and `self_induce_velocity!`.

`induce_velocity(target, source, t)` computes the complex velocity that a vortex element(s) source induces on a target at time `t`.
The target can be

- a complex position
  ```jldoctest quickstart
  julia> induce_velocity(0.0 + 0.0im , vortices, 0.0)
  0.05610938572529216 - 0.1319030126670981im

  julia> induce_velocity(0.0 + 0.0im , sys, 0.0)
  0.14592914759546077 - 0.1264803675281937im

  ```
- a vortex element
  ```jldoctest quickstart
  julia> induce_velocity(p, sys, 0.0)
  -0.004302294537820467 - 0.07805396403126988im

  ```
- an array/tuple of vortex elements
  ```jldoctest quickstart
  julia> induce_velocity(vortices, sources, 0.0)
  5-element Array{Complex{Float64},1}:
    0.06454384396015585 + 0.007898382618214123im
   0.053907048316969616 + 0.02790291832733651im
     0.0706678480701265 - 0.02711822881988212im
    0.06764122439359754 + 0.011120627964923711im
    0.07894704527850091 + 0.01178636990942516im

  julia> induce_velocity(sources, sys, 0.0)
  5-element Array{Complex{Float64},1}:
      0.1406920003437716 - 0.09680661285216872im
   -0.003388443511578354 - 0.004829334922251807im
     0.03508222004503095 - 0.10591880362792691im
      0.1221234594021031 - 0.04477695050245835im
    -0.02942885432058654 - 0.03924892616438964im

  ```

The in-place version, `induce_velocity!(velocities, targets, source, t)`, computes the velocity and writes the results into a pre-allocated data structure.
For example:
```jldoctest quickstart
julia> vel_vortices = zeros(ComplexF64, length(vortices))
5-element Array{Complex{Float64},1}:
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> induce_velocity!(vel_vortices, vortices, sources, 0.0);

julia> vel_vortices
5-element Array{Complex{Float64},1}:
  0.06454384396015585 + 0.007898382618214123im
 0.053907048316969616 + 0.02790291832733651im
   0.0706678480701265 - 0.02711822881988212im
  0.06764122439359754 + 0.011120627964923711im
  0.07894704527850091 + 0.01178636990942516im

```
To make it easier to allocate velocities for more complex collections of vortex elements, the library provides the `allocate_velocity` function:
```jldoctest quickstart
julia> vels = allocate_velocity(sys);

julia> typeof(vels)
Tuple{Array{Complex{Float64},1},Array{Complex{Float64},1}}
```
The code above created a tuple containing two arrays of velocities, corresponding to the structure of `sys`.
Similarly, there is also the `reset_velocity!(velocities, sources)` function, which resizes the entries in `velocities` to match the structure of `sources` if necessary, then sets all velocities to zero.
We can compute the velocity that a source induces on the entire points/blobs system with:
```jldoctest quickstart
julia> src = Vortex.Point(1.0, 1.0);

```

If we want the velocity that the points/blobs system induces on itself, we can call
```julia
reset_velocity!(vels, sys)
induce_velocity!(vels[1], vortices, vortices)
induce_velocity!(vels[1], vortices, sources)
induce_velocity!(vels[2], sources, vortices)
induce_velocity!(vels[2], sources, sources)
```
This becomes difficult to keep track of when `sys` gets larger or more complicated (e.g. nested collection of elements).
Instead, we can use the `self_induce_velocity!` function, which takes care of applying all the pairwise interactions (recursively if need be):
```jldoctest quickstart
julia> reset_velocity!(vels, sys);

julia> self_induce_velocity!(vels, sys, 0.0);
```

## Time Marching

```@setup timemarching
using PotentialFlow
using Plots
using Random
clibrary(:colorbrewer)
Random.seed!(1)
default(colorbar_title=("Γ"), grid = false, ratio = 1, legend = :none, colorbar = :right, markerstrokealpha = 0, markersize = 5, size = (600, 400))
```
Now that we compute the velocities of a system of vortex elements, we can march the system forward in time to simulate its behavior.
As an example, we will simulate of two clusters of vortex blobs merging.
```@example timemarching
N = 200
zs = Complex.(0.5randn(N), 0.5randn(N))
Γs  = @. exp(-4abs2(zs))
cluster₁ = Vortex.Blob.(zs .+ 1, Γs, 0.01)
cluster₂ = Vortex.Blob.(zs .- 1, Γs, 0.01)

sys = (cluster₁, cluster₂)
vels = allocate_velocity(sys)
plot(sys, color = :reds, clim = (0, 1))
savefig("initial_clusters.svg"); nothing # hide
```
```@raw html
<object data="manual/initial_clusters.svg" type="image/svg+xml"></object>
```

Given an array or tuple of vortex elements and their velocities, we can compute their positions after some time interval with the `advect!(x₊, x, ẋ, Δt)` function, where
- `x₊` is where the new states are stored
- `x` is the current state
- `Δt` is the time interval
- `ẋ` is the velocity.
In our case, we will let `x₊` and `x` both be set to `sys`:
```@example timemarching
Δt = 0.01
for t in 0:Δt:1.0
    reset_velocity!(vels, sys)
    self_induce_velocity!(vels, sys, t)
    advect!(sys, sys, vels, Δt)
end
plot(sys, color = :reds, clim = (0, 1))
savefig("final_clusters.svg"); nothing # hide
```
```@raw html
<object data="manual/final_clusters.svg" type="image/svg+xml"></object>
```
