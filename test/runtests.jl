using Test

using PotentialFlow
using LinearAlgebra: norm, ⋅

@test isempty(detect_ambiguities(Vortex))

#=
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
=#

const GROUP = get(ENV, "GROUP", "All")


if GROUP == "All" || GROUP == "Autodiff"
  include("forwarddiff.jl")
end

if GROUP == "All" || GROUP == "Bodies"
  include("Bodies.jl")
  include("utils/circle_plane.jl")
end

if GROUP == "All" || GROUP == "Plates"
  include("chebyshev.jl")
  include("Plates.jl")
end

if GROUP == "All" || GROUP == "Sheets"
  include("Sheets.jl")
end

if GROUP == "All" || GROUP == "TimeMarching"
  include("TimeMarching.jl")
end

if GROUP == "All" || GROUP == "Dispatch"
  include("dispatch.jl")
end

if GROUP == "All" || GROUP == "Utils"
  include("Utils.jl")
end

if GROUP == "All" || GROUP == "Notebooks"
  include("notebooks.jl")
end
