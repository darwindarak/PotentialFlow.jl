module Properties

using MacroTools: @capture, postwalk, striplines

export @property

## Macro to generate properties

macro property(rawexpr)
    ex = (rawexpr)

    signature = Nullable()
    stype = Nullable()
    reduce_op = Nullable()
    preallocator = Nullable()
    has_sources = has_targets = false

    postwalk(ex) do x
        @capture(x, signature = f_) && (signature = Nullable(f))
        @capture(x, stype = t_)             && (stype = Nullable(t))
        @capture(x, reduce = op_)           && (reduce_op = Nullable(op))
        @capture(x, preallocator = p_)      && (preallocator = Nullable(p))

        has_sources |= @capture(x, _::Source)
        has_targets |= @capture(x, _::Target)

        return x
    end

    if isnull(signature)
        throw(ArgumentError("missing `signature = ...` entry"))
    end

    if !(@capture(get(signature), _(__)))
        throw(ArgumentError("signature must look like a function call"))
    end

    if has_sources && has_targets
        return esc(induced_property(get(signature), stype, preallocator))
    elseif has_sources
        return esc(source_property(get(signature), stype, reduce_op))
    else
        throw(ArgumentError("missing at least one `::Source` argument in signature"))
    end
end

function source_property(signature, stype, reduce_op)
    source_ind = findfirst(x -> isa(x, Expr) && (x.args[2] == :Source), signature.args[2:end])
    source_name = signature.args[source_ind+1].args[1]

    @capture(postwalk(signature) do ex
             @capture(ex, s_::Source) && (return s)
             ex
             end, fname_(fargs__))

    @capture(postwalk(signature) do ex
             @capture(ex, s_::Source) && (return :(Elements.unwrap_targ($s)))
             ex
             end, _(unwrappedargs__))

    if isnull(reduce_op)
        mappedvars = Symbol[]
        mappedsyms = Symbol[]

        @capture(postwalk(signature) do ex
                 if @capture(ex, s_::Source)
                     sym = gensym(s)
                     push!(mappedvars, s)
                     push!(mappedsyms, sym)
                     return sym
                 end
                 ex
             end, _(mappedargs__))

        fgroup = quote
            function $fname($(fargs...), ::Type{Elements.Group})
                map($(mappedvars...)) do $(mappedsyms...)
                    $fname($(mappedargs...))
                end
            end
        end
    else
        op = get(reduce_op)

        sumvar = gensym(:Σ)
        index = gensym(:i)

        iterated_call = postwalk(signature) do ex
            @capture(ex, s_::Source) && (return :($s[$index]))
            ex
        end

        fgroup = quote
            function $fname($(fargs...), ::Type{Elements.Group})
                $sumvar = zero($(get(stype)))
                for $index in eachindex($source_name)
                    $sumvar = $op($sumvar, $iterated_call)
                end
                $sumvar
            end
        end
    end

    quote
        Core.@__doc__ function $fname($(fargs...))
            $fname($(unwrappedargs...), Elements.kind(Elements.unwrap_src($source_name)))
        end

        $fgroup
    end
end

function induced_property(signature, stype, preallocator)

    Σ = gensym(:Σ)
    index = gensym(:i)

    @capture(postwalk(signature) do ex
             @capture(ex, t_::Target) && (return t)
             @capture(ex, s_::Source) && (return s)
             ex
             end, fname_(fargs__))
    fname! = Symbol(fname, "!")

    f_allocate = get(preallocator)

    @capture(postwalk(signature) do ex
             @capture(ex, t_::Target) && (return :(Elements.unwrap_targ($t)))
             @capture(ex, s_::Source) && (return :(Elements.unwrap_targ($s)))
             ex
             end, _(unwrappedargs__))

    @capture(postwalk(signature) do ex
             @capture(ex, s_::Source) && (return :($s[$index]))
             @capture(ex, t_::Target) && (return t)
             ex
             end, _(iteratedsources__))

    @capture(postwalk(signature) do ex
             @capture(ex, s_::Source) && (return s)
             @capture(ex, t_::Target) && (return :($t[$index]))
             ex
             end, _(iteratedtargets__))

    @capture(postwalk(signature) do ex
             @capture(ex, t_::Target) && (return :($t::Tuple))
             @capture(ex, s_::Source) && (return s)
             ex
             end, _(targettuples__))


    target_ind = findfirst(x -> isa(x, Expr) && (x.args[2] == :Target), signature.args[2:end])
    target_name = signature.args[target_ind+1].args[1]

    source_ind = findfirst(x -> isa(x, Expr) && (x.args[2] == :Source), signature.args[2:end])
    source_name = signature.args[source_ind+1].args[1]

    position_args = copy(fargs)
    position_args[target_ind] = :(Elements.position($target_name))

    quote
        $f_allocate(group::Tuple) = map($f_allocate, group)

        function $f_allocate($target_name)
            $f_allocate(Elements.unwrap_targ($target_name), Elements.kind(eltype(Elements.unwrap_targ($target_name))))
        end

        $f_allocate($target_name, el::Type{Elements.Singleton}) = zeros($(get(stype)), size($target_name))

        function $fname($(fargs...))
            $fname($(unwrappedargs...),
                   Elements.kind(Elements.unwrap_targ($target_name)),
                   Elements.kind(Elements.unwrap_src($source_name)))
        end

        function $fname!(out, $(fargs...))
            $fname!(out, $(unwrappedargs...),
                    Elements.kind(Elements.unwrap_targ($target_name)),
                    Elements.kind(Elements.unwrap_src($source_name)))
        end

        function $fname($(fargs...), ::Type{Elements.Singleton}, ::Type{Elements.Singleton})
            $fname($(position_args...))
        end

        function $fname($(fargs...), ::Type{Elements.Singleton}, ::Type{Elements.Group})
            $Σ = zero($(get(stype)))
            for $index in eachindex($source_name)
                $Σ += $fname($(iteratedsources...))
            end
            $Σ
        end

        function $fname($(fargs...), ::Type{Elements.Group}, ::Any)
            out = $f_allocate($target_name)
            $fname!(out, $(fargs...))
        end

        function $fname!(out, $(fargs...), ::Type{Elements.Group}, ::Any)
            for $index in eachindex($target_name)
                out[$index] += $fname($(iteratedtargets...))
            end
            out
        end

        function $fname!(out::Tuple, $(targettuples...))
            for $index in eachindex($target_name)
                $fname!(out[$index], $(iteratedtargets...))
            end
            out
        end
    end
end

end
