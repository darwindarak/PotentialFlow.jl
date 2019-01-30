import PotentialFlow.Utils: @get, MappedVector

@testset "Utils" begin
    @testset "get" begin
        z = rand(ComplexF64)

        @get z (re, im) (r, i)
        @test i == imag(z)
        @test r == real(z)

    end

    @testset "MappedVector" begin
        x = [π, 0.0, π]
        y = MappedVector(cos, x, 1)
        @test y[0] == -1.0
        @test y[1] == 1.0
        @test y[2] == -1.0

        buff = IOBuffer()
        show(buff, y)
        @test String(take!(buff)) == "Array{Float64,1} → typeof(cos) (0:2)"
    end
end
