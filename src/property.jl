using MacroTools: @capture, postwalk, striplines
import Core: @__doc__

# begin
#     signature = induce_velocity(targ::Target, src::Source, t)
#     preallocator = allocate_velocity
#     stype = Complex128
# end

# begin
#     signature = circulation(src::Source)
#     stype = Complex128
#     reduce = +
# end

# begin
#     signature = position(src::Source)
#     stype = Complex128
# end

function source_property(signature, stype, reduce_op)
    if isnull(reduce_op)
        return :(@__doc__ function $(signature.args[1]) end)
    end

    op = get(reduce_op)

    sumvar = gensym(:Σ)
    index = gensym(:i)

    iterated_call = postwalk(signature) do ex
        @capture(ex, s_::Source) && (return :($s[$index]))
        @capture(ex, t_::Target) && (return t)
        ex
    end

    @capture(postwalk(signature) do ex
        @capture(ex, s_::Source) && (return :(unwrap_targ($s)))
        @capture(ex, t_::Target) && (return t)
        ex
    end, _(unwrappedargs__))

    @capture(postwalk(signature) do ex
        @capture(ex, s_::Source) && (return s)
        @capture(ex, t_::Target) && (return t)
        ex
    end, fname_(fargs__))
    _fname = gensym(fname)

    source_ind = findfirst(x -> x.args[2] == :Source, signature.args[2:end])
    source_name = signature.args[source_ind+1].args[1]

    quote
        @__doc__ function $fname($(fargs...))
            $_fname($(unwrappedargs...), kind(unwrap_src($source_name)))
        end

        function $_fname($(fargs...), ::Type{Group})
            $sumvar = zero($(get(stype)))
            for $index in eachindex($source_name)
                $sumvar = $op($sumvar, $iterated_call)
            end
            $sumvar
        end
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
    _fname = gensym(fname)
    _fname! = gensym(fname!)

    f_allocate = get(preallocator)
    _f_allocate = gensym(f_allocate)

    @capture(postwalk(signature) do ex
             @capture(ex, t_::Target) && (return :(unwrap_targ($t)))
             @capture(ex, s_::Source) && (return :(unwrap_targ($s)))
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


    target_ind = findfirst(x -> x.args[2] == :Target, signature.args[2:end])
    target_name = signature.args[target_ind+1].args[1]

    source_ind = findfirst(x -> x.args[2] == :Source, signature.args[2:end])
    source_name = signature.args[source_ind+1].args[1]

    position_args = copy(fargs)
    position_args[target_ind] = :(Vortex.position($target_name))

    quote
        $f_allocate(group::Tuple) = map($f_allocate, group)

        function $f_allocate($target_name)
            $_f_allocate(unwrap_targ($target_name), kind(eltype(unwrap_targ($target_name))))
        end

        $_f_allocate($target_name, el::Type{Singleton}) = zeros($(get(stype)), size($target_name))

        function $fname($(fargs...))
            $_fname($(unwrappedargs...),
                    kind(unwrap_targ($target_name)),
                    kind(unwrap_src($source_name)))
        end

        function $fname!(out, $(fargs...))
            $_fname!(out, $(unwrappedargs...),
                     kind(unwrap_targ($target_name)),
                     kind(unwrap_src($source_name)))
        end

        function $_fname($(fargs...), ::Type{Singleton}, ::Type{Singleton})
            $fname($(position_args...))
        end

        function $_fname($(fargs...), ::Type{Singleton}, ::Type{Group})
            $Σ = zero($(get(stype)))
            for $index in eachindex($source_name)
                $Σ += $fname($(iteratedsources...))
            end
            $Σ
        end

        function $_fname($(fargs...), ::Type{Group}, ::Any)
            out = $f_allocate($target_name)
            $fname!(out, $(fargs...))
        end

        function $_fname!(out, $(fargs...), ::Type{Group}, ::Any)
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

macro property(rawexpr)
    ex = striplines(rawexpr)

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
