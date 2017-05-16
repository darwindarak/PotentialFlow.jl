@testset "Vortex Sheets" begin
    @testset "Sheet construction" begin
        N = 100
        zs = rand(Complex128, N)
        Γs = accumulate(+, rand(N))
        δ = rand()

        sheet = Vortex.Sheet(zs, Γs, δ)

        @test Vortex.circulation(sheet) ≈ Vortex.circulation(sheet.blobs)
    end

    @testset "Sheet induced velocities" begin
        N = 100
        zs = rand(Complex128, N)
        Γs = accumulate(+, rand(N))
        δ = rand()
        sheet = Vortex.Sheet(zs, Γs, δ)

        ws = allocate_velocity(sheet)
        wb = allocate_velocity(sheet.blobs)

        self_induce_velocity!(ws, sheet)
        for (i, t) in enumerate(sheet.blobs), s in sheet.blobs
            wb[i] += induce_velocity(t, s)
        end
        @test wb ≈ ws
        @test wb ≈ induce_velocity(sheet, sheet)

        sheet₊ = Vortex.Sheet(Vortex.Blob[], Float64[], 0.0)
        Δt = 1e-2
        Vortex.advect!(sheet₊, sheet, ws, Δt)
        Vortex.advect!(sheet, sheet, ws, Δt)

        @test sheet₊.Γs ==  sheet.Γs
        @test !(sheet₊.Γs === sheet.Γs)

        @test sheet₊.blobs == sheet.blobs
        @test !(sheet₊.blobs === sheet.blobs)

        @test sheet₊.δ  ==  sheet.δ
    end

    @testset "Sheet surgery" begin
        N = 100
        zs = rand(Complex128, N)
        Γs = accumulate(+, rand(N))
        δ = 0.05
        sheet = Vortex.Sheet(zs, Γs, δ)

        truncate_n = rand(1:50)
        ΣΓ = Vortex.circulation(sheet)
        ΔΓ₀ = Γs[truncate_n] - Γs[1]
        ΔΓ  = Vortex.Sheets.truncate!(sheet, truncate_n)
        @test ΔΓ₀ ≈ ΔΓ

        @test Vortex.Sheets.compute_trapezoidal_weights(sheet.Γs) ≈ getfield.(sheet.blobs, :Γ)
        @test ΔΓ + Vortex.circulation(sheet) ≈ ΣΓ

        Γ = Vortex.circulation(sheet)
        sheet₊ = Vortex.Sheets.split!(sheet, rand(2:length(sheet)-2))
        @test Γ ≈ sum(Vortex.circulation, (sheet, sheet₊))

        @test Vortex.Sheets.compute_trapezoidal_weights(sheet.Γs) ≈ getfield.(sheet.blobs, :Γ)
        @test Vortex.Sheets.compute_trapezoidal_weights(sheet₊.Γs) ≈ getfield.(sheet₊.blobs, :Γ)

        N = 200
        θ = linspace(π, 0, N)
        zs = complex.(cos.(θ))
        Γs = sin.(θ)
        Vortex.Sheets.redistribute_points!(sheet, zs, Γs)
        @test Vortex.Sheets.compute_trapezoidal_weights(sheet.Γs) ≈ getfield.(sheet.blobs, :Γ)
        @test Vortex.position.(sheet.blobs) == getfield.(sheet.blobs, :z)
        @test Vortex.Sheets.positions(sheet) == Vortex.position.(sheet.blobs)

        Vortex.Sheets.remesh!(sheet, 0.005)
        @test length(sheet) == 400
        @test Vortex.Sheets.positions(sheet) ≈ linspace(-1, 1, 400)
        @test norm(sheet.Γs .- sqrt.(1 - linspace(-1,1,400).^2)) ≤ 1e-3
    end
end
