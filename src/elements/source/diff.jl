function seed(v::Vector{<:Point},cfg::ComplexGradientConfig)
  posduals = convert.(eltype(cfg.duals),Elements.position(v))
  circduals = convert.(eltype(cfg.duals),complex(flux.(v)))
  Point.(posduals,real(circduals))
end

function seed_position(v::Vector{<:Point},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v),cfg.rseeds,cfg.iseeds)
    circduals = copy(cfg.duals)
    seed!(circduals,complex(flux.(v)))
    Point.(posduals,real(circduals))
end

function seed_strength(v::Vector{<:Point},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v))
    circduals = copy(cfg.duals)
    seed!(circduals,complex(flux.(v)),cfg.rseeds,cfg.iseeds)
    Point.(posduals,real(circduals))
end

function seed(v::Vector{<:Blob},cfg::ComplexGradientConfig)
    posduals = convert.(eltype(cfg.duals),Elements.position(v))
    circduals = convert.(eltype(cfg.duals),complex(flux.(v)))
    Blob.(posduals,real(circduals),Elements.blobradius(v))
end


function seed_position(v::Vector{<:Blob},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v),cfg.rseeds,cfg.iseeds)
    circduals = copy(cfg.duals)
    seed!(circduals,complex(flux.(v)))
    Blob.(posduals,real(circduals),Elements.blobradius(v))
end

function seed_strength(v::Vector{<:Blob},cfg::ComplexGradientConfig)
    posduals = copy(cfg.duals)
    seed!(posduals,Elements.position(v))
    circduals = copy(cfg.duals)
    seed!(circduals,complex(flux.(v)),cfg.rseeds,cfg.iseeds)
    Blob.(posduals,real(circduals),Elements.blobradius(v))
end
