using Base.Test
using TestSetExtensions

using VortexModel

@test isempty(detect_ambiguities(Vortex))

@testset ExtendedTestSet "All tests" begin
    @includetests ARGS
end

if isempty(ARGS)
    include("../docs/make.jl")
end
