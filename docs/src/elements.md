# Vortex Elements

```@meta
DocTestSetup = quote
using VortexModel
srand(1)
end
```
The library currently has four built-in vortex types:

- [`Vortex.Point`](@ref)
- [`Vortex.Blob`](@ref)
- [`Vortex.Sheet`](@ref)
- [`Vortex.Plate`](@ref) (at the moment, there can only be one plate in the fluid at at time)

Most functions in the library that act on vortex elements can take either a single vortex element, or a collection of elements.
These collections can be represented as an array or a tuple.
Arrays should be used when the elements are the same type, for example:
```jldoctest overview
julia> points = Vortex.Point.(rand(Complex128, 5), rand(5))
5-element Array{VortexModel.Vortex.Points.Point,1}:
 Point Vortex: z = 0.236 + 0.347im, Γ = 0.556
 Point Vortex: z = 0.313 + 0.008im, Γ = 0.437
 Point Vortex: z = 0.489 + 0.211im, Γ = 0.425
 Point Vortex: z = 0.952 + 1.0im, Γ = 0.773
 Point Vortex: z = 0.252 + 0.987im, Γ = 0.281

julia> Vortex.impulse(points)
1.3362266530178137 - 1.2821936908564113im

julia> blobs = [Vortex.Blob(rand(Complex128), rand(), 0.1) for i in 1:5]
5-element Array{VortexModel.Vortex.Blobs.Blob,1}:
 Vortex Blob: z = 0.209 + 0.251im, Γ = 0.02, δ = 0.1
 Vortex Blob: z = 0.288 + 0.86im, Γ = 0.077, δ = 0.1
 Vortex Blob: z = 0.64 + 0.874im, Γ = 0.279, δ = 0.1
 Vortex Blob: z = 0.751 + 0.645im, Γ = 0.078, δ = 0.1
 Vortex Blob: z = 0.848 + 0.086im, Γ = 0.553, δ = 0.1

julia> Vortex.impulse(blobs)
0.41217890550975256 - 0.7325028967929701im
```
Knowing that every element has the same type allows the compiler to perform more aggressive optimizations.
Tuples are used when we want to mix and match *different* vortex types.
For example:
```julia
julia> sys = (points, blobs);

julia> Vortex.impulse(sys)
1.7484055585275664 - 2.0146965876493814im
```

This rest of this page documents the data types that represent these elements and some key functions that act on them.
For more detailed examples, please refer to the [Jupyter notebooks](https://github.com/darwindarak/VortexModel.jl/tree/master/examples).


## Built-in Vortex Types

```@docs
Vortex.Point
Vortex.Blob
Vortex.Sheet
Vortex.Plate
```

## Vortex Properties

```@docs
Vortex.position
Vortex.circulation
Vortex.impulse
advect
advect!
```

## Methods on Vortex Sheets

```@docs
Vortex.Sheets.append_segment!
Vortex.Sheets.truncate!
Vortex.Sheets.redistribute_points!
Vortex.Sheets.remesh
Vortex.Sheets.remesh!
Vortex.Sheets.split!
Vortex.Sheets.filter!
Vortex.Sheets.filter_position!
Vortex.Sheets.arclength
Vortex.Sheets.arclengths
```

## Methods on Plates

```@docs
Vortex.Plates.enforce_no_flow_through!
Vortex.Plates.vorticity_flux
Vortex.Plates.vorticity_flux!
Vortex.Plates.bound_circulation
Vortex.Plates.bound_circulation!
```

## Index

```@index
Pages = ["elements.md"]
```
