
using ForwardDiff
using DiffRules

import ForwardDiff: Dual, Partials, single_seed, partials, Chunk, Tag, seed!, value,
              gradient, valtype, extract_gradient!, derivative,extract_derivative,
              checktag, vector_mode_gradient, vector_mode_jacobian,
              chunk_mode_gradient

const AMBIGUOUS_TYPES = (AbstractFloat, Irrational, Integer, Rational,
                         Real, RoundingMode, ComplexF64)


include("dual.jl")
include("config.jl")
include("derivative.jl")
include("gradient.jl")
include("jacobian.jl")
include("apiutils.jl")
include("diffrules.jl")
