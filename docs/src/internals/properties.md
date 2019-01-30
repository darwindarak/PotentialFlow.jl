# Handing Pairwise Interactions

```@meta
DocTestSetup = quote
	using PotentialFlow
end
```

We want users to be able to define their own vortex types, as well as arbitrarily group and nest different vortex elements together.
For example, suppose the user has defined a new element, `MyVortexType`, then they should be able to do something like
```julia
myvortices = MyVortexType.(...)
points = Vortex.Point.(...)
sheet = Vortex.Sheet(...)

system = (myvortices, (points, sheet))
velocities = allocate_velocity(system)

self_induce_velocity!(velocities, system)
```
But how would `myvortices` know how to induce velocities on `points`, `sheet`, or the tuple `(points, sheet)`?
It would be asking a lot for the user to have to define all possible pairwise interactions between their new vortex type and all other built-in types.
Instead, the user should only have to define `induce_velocity` between `MyVortexType` and a complex point, leaving library to simply just apply `induce_velocity` recursively to all ther targets.
But how would the library know if a vortex element is actually a collection of more primitive types that should be recursed over?
For example, while it is obvious that `Vector{Vortex.Point}` should be treated as a collection of point vortices, it has no prior knowledge of how `MyVortexType` is defined.
It might be something like `Vortex.Blob`, which cannot recursed into, or it can be more like `Vortex.Sheet`, which is just a wrapper around `Vector{Vortex.Blob}`.
The following section describes how the library handles this problem using [Tim Holy's trait trick](https://github.com/JuliaLang/julia/issues/2345#issuecomment-54537633).

## Traits: `Singleton` vs. `Group`

Let's trace through how the library currently handles a call like
```julia
induce_velocity(target::V1, source::V2)
```
If the user has explicitly defined `induce_velocity` between the vortex types `V1` and `V2`, then Julia will call that method.
Otherwise, this call will be turned into
```julia
induce_velocity(unwrap_targ(target), unwrap_src(source),
                kind(unwrap_targ(target)), kind(unwrap_src(source)))
```
There are two different things going on here:
- `unwrap_targ` and `unwrap_src`: There are some vortex types that are essentually wrappers around more primitive types.
  For example, `Vortex.Sheet` is a wrapper around `Vector{Vortex.Blob}`, with some extra information to maintain connectivity between the blobs.
  Instead of having to redefine all the required functions for `Vortex.Sheet`, we simply define the function `Vortex.unwrap(s::Sheet) = s.blobs`.
  Then whenever the library encounters a `Vortex.Sheet`, it will know to unwrap it and treat it as an array of `Vortex.Blob`.
  By default, `unwrap_targ(v) = v` and `unwrap_src(v) = v`.

- `kind`: This is a trait function takes a vortex element and returns either `Type{Singleton}` or `Type{Group}`.
  A vortex with trait `Type{Singleton}` tells the library that it should be treated as a single entity, and should not be recursed into.
  Alternatively, an element with `Type{Group}` trait tells that library that this element is indexable and should be iterated through.

There are four possible `(target,source)` trait combinations:
```julia
induce_velocity(uw_target, uw_source, ::Type{Singleton}, ::Type{Singleton})
induce_velocity(uw_target, uw_source, ::Type{Group}, ::Type{Singleton})
induce_velocity(uw_target, uw_source, ::Type{Singleton}, ::Type{Group})
induce_velocity(uw_target, uw_source, ::Type{Group}, ::Type{Group})
```
but we only have to handle three cases:
- `(::Type{Singleton}, ::Type{Singleton})`:
  The fact that the call chain got to this point at all means, that there is no specialized `induce_velocity` defined between `uw_target` and `uw_source`, otherwise Julia's dispatch system would have call that one instead (see the `induce_velocity` definitions in `Plates.jl`).
  Since all vortex types are required to define `induced_velocity` on a point, this call is turned into
  ```julia
  induce_velocity(Vortex.position(uw_target), uw_soruce)
  ```

- `(::Type{Singleton}, ::Type{Group})`: In this case, we iterate through `i in 1:length(uw_source)` and sum up the the results of `induce_velocity(uw_target, uw_source[i])`

- `(::Type{Group}, ::Any)`: Since the output is no longer a scalar value, we first preallocate the output with `allocate_velocity(uw_target)`, then iteratively apply `induce_velocity` over all the elements of `uw_target`.
  Once the target has been expanded all the way to `Singleton` types, then we are back to the `(target, source)` kind being either `(::Type{Singleton}, ::Type{Group})` or `(::Type{Singleton}, ::Type{Singleton})`, which can be handled by the two cases listed above.

Ultimately, this whole setup is just a way to allow `induce_velocity` to be called recursively on arbitrary groupings of vortex elements.
However, velocity is not the only property that can be computed this way.
Other quantities, such as acceleration, circulation, etcs., can all be be computed using the same framework.
Instead of writing essentially the same code for all of them, we can use the `@property` macro

## The `@property` macro

All the `induce_velocity` methods listed above (and their in-place version, `induce_velocity!`) can be generated with
```julia
@property begin
    signature = induce_velocity(targ::Target, src::Source)
    preallocator = allocate_velocity
    stype = ComplexF64
end
```
where
- `signature` tells the macro what the function should be called as well as the roles of the different arguments.
- `preallocator` is the name you want to use to allocate output space for the targets
- `stype` is the data type of the output.

!!! note
    You can see the actual functions generated with
    ```julia
    julia> import VortexModel.Vortex: @property

    julia> @macroexpand(@property begin
               signature = induce_velocity(targ::Target, src::Source)
               preallocator = allocate_velocity
               stype = ComplexF64
           end)
    ```

Suppose we want a function to also get the acceleration of the vortex elements.
To find the acceleration, we need to know the current velocity, so we can have something like
```julia
@property begin
    signature = induce_acceleration(targ::Target, targvel::Target, src::Source, srcvel::Source)
    preallocator = allocate_acceleration
    stype = ComplexF64
end
```
By annotating `targvel` as a `Target`, we are saying that whenever you iterate through `targ`, we should pass in the corresponding element in `targvel`, and likewise for `srcvel`.
Arguments that are not annotated will be treated as parameters to be passed in at each iteration without indexing.

There are also properties that do not require a target, but are properties of the source itself.
For example, we have
```julia
@property begin
    signature = circulation(src::Source)
    stype = Float64
    reduce = (+)
end
```
The `reduce` operation means that it is a property that can be aggregate over a collection of vortex elements.
In this particular case, it means that the circulation of a group of vortex elements is just the sum of the circulation of each element.
Another example of this can be seen in the definition of [`rate_of_impulse`](@ref Plates.rate_of_impulse).

### Defining a new property

We'll go through an example of how to define new properties using the `@property` marco.
Suppose we want to check if a system of elements have branch cuts in their streamfunction, we can simply define the following:
```jldoctest
import PotentialFlow.Properties: @property

@property begin
	signature = continuous_streamfunction(src::Source)
	stype = Bool
	reduce = (&, true)
end

continuous_streamfunction(::Vortex.Point) = true
continuous_streamfunction(::Vortex.Blob) = true
continuous_streamfunction(::Source.Point) = false
continuous_streamfunction(::Source.Blob) = false

vortices = (Vortex.Point.(rand(10), rand(10)), 
	        Vortex.Blob.(rand(10), rand(10), rand())
		   )

sources = (Source.Point.(rand(10), rand(10)),
           Source.Blob.(rand(10), rand(10), rand())
		  )

mixed = (vortices, sources)

continuous_streamfunction.((vortices, sources, mixed))

# output

(true, false, false)
```
Here, the `reduce` operation is a tuple that takes in a binary operation and an initial value.
When `continuous_streamfunction` is called on a group source, such as an array of elements, it will recursively call `continuous_streamfunction` on each member of the group, and use `&` to combine the results.
Without the `true` initial value, the `@property` macro will use `zero(stype)`, which in this case, would have been `false`.
If we did not want the values to be aggregated, but instead wanted to preserve the organization structure of our source elements, we can simply leave out the `reduce` field.
For instance, if we wanted to know whether the element is a desingularized element or not, it does not make sense to reduce the results.
```jldoctest
import PotentialFlow.Properties: @property
import PotentialFlow: Points, Blobs

@property begin
	signature = is_desingularized(src::Source)
	stype = Bool
end

is_desingularized(::Points.Point) = false
is_desingularized(::Blobs.Blob) = true

vortices = (Vortex.Point.(rand(2), rand(2)), 
	        Vortex.Blob.(rand(2), rand(2), rand())
		   )

sources = (Source.Point.(rand(2), rand(2)),
           Source.Blob.(rand(2), rand(2), rand())
		  )

is_desingularized.((vortices, sources))

# output

((Bool[false, false], Bool[true, true]), (Bool[false, false], Bool[true, true]))
```
