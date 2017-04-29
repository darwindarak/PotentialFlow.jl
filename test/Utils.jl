@testset "Utils" begin
    @testset "get" begin
        z = rand(Complex128)

        @get z (im, re)
        @test im == imag(z)
        @test re == real(z)

        @get z (re, im)
        @test im == imag(z)
        @test re == real(z)

        @get z (re, im) (r, i)
        @test i == imag(z)
        @test r == real(z)

        @test_throws ErrorException (@eval @get z z)
        @test_throws ErrorException (@eval @get z (re, im) nothing)
        @test_throws ErrorException (@eval @get z (re, im) (r, i) nothing)
        @test_throws AssertionError (@eval @get z (re, im) (im,))
        @test_throws ErrorException (@eval @get z (re, im) (im,) (im,))
    end

    @testset "MappedVector" begin
        x = [π, 0.0, π]
        y = Vortex.Utils.MappedVector(cos, x, 1)
        @test y[0] == -1.0
        @test y[1] == 1.0
        @test y[2] == -1.0

        buff = IOBuffer()
        show(buff, y)
        @test String(take!(buff)) == "Array{Float64,1} → Base.#cos (0:2)"
    end
end
