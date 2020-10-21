
#=
seed is designed to form duals of the positions and strengths of the
type given in cfg
=#
function seed(v::Vector{<:Point},cfg::ComplexGradientConfig)
  posduals = convert.(eltype(cfg.duals),Elements.position(v))
  circduals = convert.(eltype(cfg.duals),complex(circulation.(v)))
  Point.(posduals,real(circduals))
end
function seed(v::Vector{<:Blob},cfg::ComplexGradientConfig)
  posduals = convert.(eltype(cfg.duals),Elements.position(v))
  circduals = convert.(eltype(cfg.duals),complex(circulation.(v)))
  Blob.(posduals,real(circduals),Elements.blobradius(v))
end

#=
seed_position will create complex duals with of the positions, with 1s in
the partials. It will also create complex duals of the strengths, even though
these are real-valued, since there is the possibility for strengths'
partials to become complex (e.g. via edge conditions), and we need to ensure
that operations on these duals get dispatched to our suite of ComplexDual tools.
=#
function seed_position(v::Vector{<:Point},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v),cfg.rseeds,cfg.iseeds)
    circduals = copy(cfg.duals)
    seed!(circduals,complex(circulation.(v)))
    Point.(posduals,real(circduals))
end
function seed_position(v::Vector{<:Blob},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v),cfg.rseeds,cfg.iseeds)
    circduals = copy(cfg.duals)
    seed!(circduals,complex(circulation.(v)))
    Blob.(posduals,real(circduals),Elements.blobradius(v))
end

#=
seed_strength will create complex duals with of the strengths, with 1s in
the partials. It does this even though the strengths are themselves real,
since there is the possibility for strengths' partials to become complex
(e.g. via edge conditions), and we need to ensure
that operations on these duals get dispatched to our suite of ComplexDual tools.
It will also create complex duals of the positions.
=#
function seed_strength(v::Vector{<:Point},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v))
    circduals = copy(cfg.duals)
    seed!(circduals,complex(circulation.(v)),cfg.rseeds,cfg.iseeds)
    Point.(posduals,real(circduals))
end
function seed_strength(v::Vector{<:Blob},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v))
    circduals = copy(cfg.duals)
    seed!(circduals,complex(circulation.(v)),cfg.rseeds,cfg.iseeds)
    Blob.(posduals,real(circduals),Elements.blobradius(v))
end
