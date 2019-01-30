using FFTW

@testset "Chebyshev Routines" begin
    import PotentialFlow.Plates.Chebyshev

    @testset "Quadrature" begin
        N = 127
        w = Chebyshev.clenshaw_curtis_weights(N)
        @test sum(w) ≈ 2

        x = Chebyshev.nodes(N)
        @test w ⋅ x ≤ eps()
        @test w ⋅ x.^2 ≈ 2/3
        @test w ⋅ exp.(x) ≈ (exp(1) - exp(-1))
    end

    @testset "Series" begin
        # Check evaluation off offset Chebyshev series
        s = rand(2)
        @test Chebyshev.firstkind([rand(), 2.0, 1.0], s, 1)     ≈ @. 2s^2 + 2s - 1
        @test Chebyshev.firstkind([rand(), rand(), 1.0], s, 2)  ≈ @. 2s^2 - 1
        @test Chebyshev.secondkind([rand(), 2.0, 1.0], s, 1)    ≈ @. 4s^2 + 4s - 1
        @test Chebyshev.secondkind([rand(), rand(), 1.0], s, 2) ≈ @. 4s^2 - 1
    end

    @testset "Transforms" begin

        N = rand(128:512)
        R = zeros(Float64, N)
        y = zeros(Float64, N)
        @test all(rand(0:(N÷2), 10)) do n
            x = [cos(n*θ) for θ in range(π, 0, length=N)]
            Chebyshev.transform!(R, x)
            Chebyshev.inv_transform!(y, R)
            (y ≈ x) && (R[n+1] ≈ 1.0) && (sum(R) ≈ 1.0)
        end

        n₁ = rand(0:(N÷2))
        n₂ = rand(0:(N÷2))
        x₁ = [cos(n₁*θ) for θ in range(π, 0, length=N)]
        x₂ = [cos(n₂*θ) for θ in range(π, 0, length=N)]

        C = x₁ .+ im.*x₂
        I = zeros(Float64, N)

        Chebyshev.transform!(R, x₁)
        Chebyshev.transform!(I, x₂)
        Chebyshev.transform!(C)

        @test real.(C) ≈ R
        @test imag.(C) ≈ I

        plan! = FFTW.plan_r2r!(C, FFTW.REDFT00)
        @allocated Chebyshev.transform!(C, plan!)
        @test 0 == (@allocated Chebyshev.transform!(C, plan!))

        @testset "Planned Transforms" begin
            x = Chebyshev.nodes(N)

            K  = Chebyshev.plan_transform(ones(N))
            K! = Chebyshev.plan_transform!(ones(N))
            @test_throws BoundsError K*ones(N+1)
            @test_throws BoundsError K\ones(N+1)

            C = K*ones(N)
            @test C[1] ≈ 1
            @test norm(C[2:end]) ≤ 128eps()

            C = K*x
            @test norm(C[1]) ≤ 128eps()
            @test C[2] ≈ 1
            @test norm(C[3:end]) ≤ 128eps()

            y = @. 16x^5 - 20x^3 + 5x
            C = K*y
            @test norm(C[1:5]) ≤ 128eps()
            @test C[6] ≈ 1
            @test norm(C[7:end]) ≤ 128eps()

            ỹ = K \ C
            @test y ≈ ỹ

            # Inverse transforms should be the same as simply
            # evaluating the series
            @test Chebyshev.firstkind(C, x) ≈ ỹ

            K! \ C
            @test y ≈ C

            # Make sure transforms work with even number of nodes
            N = 128
            x = Chebyshev.nodes(N)
            K! = Chebyshev.plan_transform!(x)
            y = @. 16x^5 - 20x^3 + 5x
            C = copy(y)
            K! * C
            ỹ = Chebyshev.inv_transform(C)
            @test y ≈ ỹ

            K! \ C
            @test y ≈ C

            # Make sure transforms work with odd number of nodes
            N = 127
            x = Chebyshev.nodes(N)
            K! = Chebyshev.plan_transform!(x)
            y = @. 16x^5 - 20x^3 + 5x
            C = copy(y)
            K! * C
            ỹ = Chebyshev.inv_transform(C)
            @test y ≈ ỹ

            K! \ C
            @test y ≈ C
        end
    end

end
