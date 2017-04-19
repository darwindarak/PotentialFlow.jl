module Utils

export @submodule, @get

macro submodule(mod)
    path = mod * ".jl"
    name = split(mod, '/')[end]
    Expr(:block,
         Expr(:call, :include, path),
         Expr(:using, :., Symbol(name)))
end

"""
A macro for extracting fields from an object.  For example, instead of a statement
like

    (obj.a + obj.b)^obj.c

it is more readable if we first locally bind the variables

    a = obj.a
    b = obj.b
    c = obj.c

    (a + b)^c

Using the `@get` marco, this becomes

    @get obj (a, b, c)

    (a + b)^c

We can also locally assign different names to the binding

    @get obj (a, b, c) (α, β, γ)
    (α + β)^γ
"""
macro get(object, fields...)
    if length(fields) == 1
        try
            @assert typeof(fields[1]) == Expr
            @assert fields[1].head == :tuple
            @assert all([typeof(arg) == Symbol for arg in fields[1].args])
        catch
            error("second argument must be a tuple of field names")
        end
        esc(Expr(:block, [:($sym = $object.$sym) for sym in fields[1].args]...))
    elseif length(fields) == 2
        try
            @assert typeof(fields[1]) == Expr
            @assert typeof(fields[2]) == Expr
            @assert (fields[1].head == :tuple && fields[2].head == :tuple)
            @assert all([typeof(arg) == Symbol for arg in fields[1].args])
            @assert all([typeof(arg) == Symbol for arg in fields[2].args])
        catch
            error("second and third argument must be tuples of field names")
        end

        nargs = length(fields[1].args)
        @assert nargs == length(fields[2].args) "field name tuples must have the same length"
        esc(Expr(:block, [:($(fields[2].args[i]) = $object.$(fields[1].args[i])) for i in 1:nargs]...))
    else
        error("Usage: @get <object> (field names...) [(reference names...)]")
    end
end

end
