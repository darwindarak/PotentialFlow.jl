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
5-element Vector{PotentialFlow.Points.Point{Float64, Float64, Val{Inf}}}:
 Vortex.Point(0.07336635446929285 + 0.34924148955718615im, 0.5710874493423871)
 Vortex.Point(0.6988266836914685 + 0.6282647403425017im, 0.4528085872833483)
 Vortex.Point(0.9149290036628314 + 0.19280811624587546im, 0.30232547191787174)
 Vortex.Point(0.7701803478856664 + 0.7805192636751863im, 0.0013502779247226426)
 Vortex.Point(0.6702639583444937 + 0.16771210647092682im, 0.5670236732404312)

julia> Elements.impulse(points)
0.638372558313404 - 1.0160351596663442im

julia> blobs = [Vortex.Blob(rand(ComplexF64), rand(), 0.1) for i in 1:5]
5-element Vector{PotentialFlow.Blobs.Blob{Float64, Float64, Val{Inf}}}:
 Vortex.Blob(0.6159379234562881 + 0.19573857852575793im, 0.012461945950411835, 0.1)
 Vortex.Blob(0.3119923865097316 + 0.11479916823306191im, 0.5460487092960259, 0.1)
 Vortex.Blob(0.6232150941621899 + 0.2708693898950604im, 0.8451820156319791, 0.1)
 Vortex.Blob(0.49359045543272007 + 0.9003405842788204im, 0.37215957409032674, 0.1)
 Vortex.Blob(0.8686942572391998 + 0.8667711192602672im, 0.7305508461555176, 0.1)

julia> Elements.impulse(blobs)
1.2623499011326258 - 1.523088752876448im
```
Knowing that every element has the same type allows the compiler to perform more aggressive optimizations.
Tuples are used when we want to mix and match *different* element types.
For example:
```julia
julia> sys = (points, blobs);

julia> Elements.impulse(sys)
1.9007224594460297 - 2.5391239125427925im
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
