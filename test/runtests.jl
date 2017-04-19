using Base.Test
using TestSetExtensions

using VortexModel

@test isempty(detect_ambiguities(Vortex))

test_files = [
    "dispatch.jl",
    "Plates.jl",
    "Sheets.jl",
]

@testset DottedTestSet "All tests" begin
    @includetests ARGS
end
