using Base.Test
using TestSetExtensions

using VortexModel

@test isempty(detect_ambiguities(Vortex))

@testset DottedTestSet "All tests" begin
    @includetests ARGS
end
