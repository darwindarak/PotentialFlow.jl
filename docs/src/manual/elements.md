# Elements

```@meta
DocTestSetup = quote
using PotentialFlow
using Random
Random.seed!(1)
end
```
The library currently has these built-in potential flow elements:

- [`Vortex.Point`](@ref)
- [`Vortex.Blob`](@ref)
- [`Vortex.Sheet`](@ref)
- [`Source.Point`](@ref)
- [`Source.Blob`](@ref)
- [`Plate`](@ref) (at the moment, there can only be one plate in the fluid at at time)
- [`Bodies.ConformalBody`](@ref)

Most functions in the library that act on elements can take either a single element, or a collection of elements.
These collections can be represented as an array or a tuple.
Arrays should be used when the elements are the same type, for example:
```jldoctest overview
julia> points = Vortex.Point.(rand(ComplexF64, 5), rand(5))
5-element Array{PotentialFlow.Points.Point{Float64},1}:
 Vortex.Point(0.23603334566204692 + 0.34651701419196046im, 0.5557510873245723)
 Vortex.Point(0.3127069683360675 + 0.00790928339056074im, 0.43710797460962514)
 Vortex.Point(0.4886128300795012 + 0.21096820215853596im, 0.42471785049513144)
 Vortex.Point(0.951916339835734 + 0.9999046588986136im, 0.773223048457377)
 Vortex.Point(0.25166218303197185 + 0.9866663668987996im, 0.2811902322857298)

julia> Elements.impulse(points)
1.3362266530178137 - 1.2821936908564113im

julia> blobs = [Vortex.Blob(rand(ComplexF64), rand(), 0.1) for i in 1:5]
5-element Array{PotentialFlow.Blobs.Blob{Float64},1}:
 Vortex.Blob(0.20947237319807077 + 0.25137920979222494im, 0.02037486871266725, 0.1)
 Vortex.Blob(0.2877015122756894 + 0.859512136087661im, 0.07695088688120899, 0.1)
 Vortex.Blob(0.6403962459899388 + 0.8735441302706854im, 0.27858242002877853, 0.1)
 Vortex.Blob(0.7513126327861701 + 0.6448833539420931im, 0.07782644396003469, 0.1)
 Vortex.Blob(0.8481854810000327 + 0.0856351682044918im, 0.5532055454580578, 0.1)

julia> Elements.impulse(blobs)
0.41217890550975256 - 0.7325028967929701im
```
Knowing that every element has the same type allows the compiler to perform more aggressive optimizations.
Tuples are used when we want to mix and match *different* element types.
For example:
```julia
julia> sys = (points, blobs);

julia> Elements.impulse(sys)
1.7484055585275664 - 2.0146965876493814im
```

This rest of this page documents the data types that represent these elements and some key functions that act on them.
For more detailed examples, please refer to the [Jupyter notebooks](https://github.com/darwindarak/PotentialFlow.jl/tree/master/examples).

## Built-in Types

```@docs
Vortex.Point
Vortex.Blob
Vortex.Sheet
Source.Point
Source.Blob
Plate
Bodies.ConformalBody
```

## Element Properties

```@docs
Elements.position
Elements.circulation
Elements.flux
Elements.impulse
```

## Methods on Vortex Sheets

```@docs
Sheets.append_segment!
Sheets.truncate!
Sheets.redistribute_points!
Sheets.remesh
Sheets.remesh!
Sheets.split!
Sheets.filter!
Sheets.filter_position!
Sheets.arclength
Sheets.arclengths
```

## Methods on Plates

```@docs
Plates.edges
Plates.enforce_no_flow_through!
Plates.vorticity_flux
Plates.vorticity_flux!
Plates.bound_circulation
Plates.bound_circulation!
Plates.rate_of_impulse
Plates.force
Plates.surface_pressure
```

## Methods on Conformally-Mapped Bodies

```@docs
Bodies.enforce_no_flow_through!
Bodies.normal
Bodies.tangent
Bodies.transform_velocity!
```

## Index

```@index
Pages = ["elements.md"]
```
