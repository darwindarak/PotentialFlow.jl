using MacroTools

"""
    @property kind name property_type

Macro to define common functions for computing vortex properties

`kind` can be one of:
- `induced`: when the property is induced by an external vortex source (e.g. velocity),
  the macro defines the following functions:
  ```julia
  allocate_<name>(target)
  induce_<name>(target, source)
  induce_<name>!(output, target, source)
  ```
  where `source` can be
  - a single vortex element (e.g. a single point vortex)
  - a array of homogenous vortex elements (e.g. a group of point vortices)
  - a tuple of different vortex elements (e.g. a plate and two vortex sheets),
  `target` can be the same types as `source` as well as one or more positions,
  and `output` is the output array allocated by `allocate_<name>`
  The return type of `induce_<name>` will be `property_type`.
- `point`: when the property is pointwise defined (e.g. position),
  the macro defines a placeholder function
  ```julia
  function <name> end
  ```
  deferring the actual implementation to the indivual vortex types.
- `aggregate`: when the property can be summed over multiple elements (e.g. circulation),
  the macro defines
  ```julia
  <name>(source)
  ```

!!! note Why a macro?
    Originally, `induce_velocity` and friends were defined as regular
    functions.  However, computations of many other properties
    (e.g. complex potential, acceleration, etc.)  follow the same
    pattern of summing over sources and distributing over targets.
    Instead of writing separate code, it seems better to write one set
    of code that will work on all properties of this kind.  That way,
    defining new properties in the future will be eaiser, and we only
    have one set of bugs to fix.
"""
macro property(kind, name, prop_type, target_deps = :(()), source_deps = :(()))
    if kind == :point
        esc(quote
            function $name end
            end)
    elseif kind == :aggregate
        esc(quote
            $name(vs::Collection) = mapreduce($name, +, zero($prop_type), vs)
            end)
    elseif kind == :induced
        f_allocate = Symbol("allocate_$(lowercase(string(name)))")
        f_induce   = Symbol("induce_$(lowercase(string(name)))")
        f_induce!  = Symbol("induce_$(lowercase(string(name)))!")

        @capture(target_deps, (t_deps__,))
        @capture(source_deps, (s_deps__,))
        t_deps_i = map(s -> :($s[i]), t_deps)
        s_deps_i = map(s -> :($s[i]), s_deps)

        esc(quote
            function $f_allocate(array::AbstractArray{T}) where {T <: Union{PointSource, Complex128}}
                zeros($prop_type, size(array))
            end

            $f_allocate(group::Tuple) = map($f_allocate, group)

            function $f_induce(target::PointSource, $(t_deps...), source, $(s_deps...))
                $f_induce(Vortex.position(target), $(t_deps...), source, $(s_deps...))
            end

            function $f_induce(z::Complex128, $(t_deps...), sources::Collection, $(s_deps...))
                w = zero($prop_type)
                for i in eachindex(sources)
                    w += $f_induce(z, $(t_deps...), sources[i], $(s_deps_i...))
                end
                w
            end

            function $f_induce(targets::Collection, $(t_deps...), source, $(s_deps...))
                ws = $f_allocate(targets)
                $f_induce!(ws, targets, $(t_deps...), source, $(s_deps...))
            end

            function $f_induce!(ws::AbstractArray, targets::AbstractArray, $(t_deps...), source, $(s_deps...))
                for i in eachindex(targets)
                    ws[i] += $f_induce(targets[i], $(t_deps_i...), source, $(s_deps...))
                end
                ws
            end

            function $f_induce!(ws::Tuple, targets::Tuple, $(t_deps...), source, $(s_deps...))
                for i in 1:length(targets)
                    $f_induce!(ws[i], targets[i], $(t_deps_i...), source, $(s_deps...))
                end
                ws
            end
            end)
    else
        throw(ArgumentError("first argument must be `:induced`, `:point` or `:aggregate`"))
    end
end
