using Test

using PotentialFlow
using LinearAlgebra: norm, ⋅ 

@test isempty(detect_ambiguities(Vortex))

# Copied verbatim from ssfrr/TestSetExtensions.jl
# The package is not 0.7/1.0 compatible yet, but
# this function is too useful
macro includetests(testarg...)
    if length(testarg) == 0
        tests = []
    elseif length(testarg) == 1
        tests = testarg[1]
    else
        error("@includetests takes zero or one argument")
    end

    quote
        tests = $tests
        rootfile = @__FILE__
        if length(tests) == 0
            tests = readdir(dirname(rootfile))
            tests = filter(f->endswith(f, ".jl") && f!= basename(rootfile), tests)
        else
            tests = map(f->string(f, ".jl"), tests)
        end
        println();
        for test in tests
            print(splitext(test)[1], ": ")
            include(test)
            println()
        end
    end
end

@testset "All tests" begin
    @includetests ARGS
end
