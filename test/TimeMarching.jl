@testset "Time Marching" begin
    import PotentialFlow.Utils: centraldiff
    T = [1e-1, 1e-2, 1e-3]

    f! = (ẋ, x, t) -> ẋ .= [-x[2], x[1]]
    update = (x₊, x₋, ẋ, Δt) -> @. x₊ = x₋ + Δt*ẋ

    @testset "Forward Euler" begin
        err = map(T) do Δt
            x₋ = [1.0, 2.0]
            x₊ = similar(x₋)
            ẋ  = similar(x₋)

            for t in 0:Δt:5.0-Δt
                forward_euler!(x₊, x₋, t, Δt, f!, update, ẋ)
                x₊, x₋ = x₋, x₊
            end
            norm(x₋ - [cos(5)-2sin(5),2cos(5)+sin(5)])
        end

        @test minimum(centraldiff(log.(err))./centraldiff(log.(T))) > 0.9
    end

    @testset "RK4" begin
        err = map(T) do Δt
            x₋ = [1.0, 2.0]
            x₊ = similar(x₋)
            ẋ  = [similar(x₋) for k in 1:4]

            for t in 0:Δt:5.0-Δt
                rk4!(x₊, x₋, t, Δt, f!, update, ẋ)
                x₊, x₋ = x₋, x₊
            end
            norm(x₋ - [cos(5)-2sin(5),2cos(5)+sin(5)])
        end

        @test minimum(centraldiff(log.(err))./centraldiff(log.(T))) > 3.9
    end
end
