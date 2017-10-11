module Elements

export Element, Singleton, Group, kind, @kind, circulation

using ..Properties

abstract type Element end

#== Trait definitions ==#

abstract type Singleton end
abstract type Group end

function kind end

macro kind(element, k)
    esc(quote
        Elements.kind(::$element) = $k
        Elements.kind(::Type{$element}) = $k
    end)
end

@kind Complex128 Singleton
kind(::AbstractArray{T}) where {T <: Union{Element, Complex128}} = Group
kind(::Tuple) = Group

# Convenience functions to define wrapper types
# e.g. vortex sheets as a wrapper around vortex blobs
unwrap_src(e) = unwrap(e)
unwrap_targ(e) = unwrap(e)
unwrap(e) = e

## Actual Definitions of Properties

@property begin
    signature = position(src::Source)
    stype = Complex128
end

@property begin
    signature = circulation(src::Source)
    reduce = (+)
    stype = Float64
end

@property begin
    signature = impulse(src::Source)
    reduce = (+)
    stype = Complex128
end

end
