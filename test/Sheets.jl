@testset "Vortex Sheets" begin
    @testset "Sheet construction" begin
        N = 100
        zs = rand(ComplexF64, N)
        Γs = accumulate(+, rand(N))
        δ = rand()

        sheet = Vortex.Sheet(zs, Γs, δ)

        @test circulation(sheet) ≈ circulation(sheet.blobs)
    end

    @testset "Sheet induced velocities" begin
        N = 100
        zs = rand(ComplexF64, N)
        Γs = accumulate(+, rand(N))
        δ = rand()
        sheet = Vortex.Sheet(zs, Γs, δ)

        ws = allocate_velocity(sheet)
        wb = allocate_velocity(sheet.blobs)

        self_induce_velocity!(ws, sheet, 0)
        for (i, t) in enumerate(sheet.blobs), s in sheet.blobs
            wb[i] += induce_velocity(t, s, 0)
        end
        @test wb ≈ ws
        @test wb ≈ induce_velocity(sheet, sheet, 0)

        sheet₊ = Vortex.Sheet(Vector{Vortex.Blob}(undef, 10), rand(10), 0.0)
        Δt = 1e-2
        advect!(sheet₊, sheet, ws, Δt)
        advect!(sheet, sheet, ws, Δt)

        @test sheet₊.Ss ==  sheet.Ss
        @test !(sheet₊.Ss === sheet.Ss)

        @test sheet₊.blobs == sheet.blobs
        @test !(sheet₊.blobs === sheet.blobs)

        @test sheet₊.δ  ==  sheet.δ
    end

    @testset "Sheet surgery" begin
        N = 100
        zs = rand(ComplexF64, N)
        Γs = accumulate(+, rand(N))
        δ = 0.05
        sheet = Vortex.Sheet(zs, Γs, δ)

        truncate_n = 50
        ΣΓ = circulation(sheet)
        ΔΓ₀ = Γs[truncate_n] - Γs[1]
        ΔΓ  = Sheets.truncate!(sheet, truncate_n)
        @test ΔΓ₀ ≈ ΔΓ

        @test Sheets.compute_trapezoidal_weights(sheet.Ss) ≈ circulation.(sheet.blobs)
        @test ΔΓ + circulation(sheet) ≈ ΣΓ

        Γ = circulation(sheet)
        sheet₊ = Sheets.split!(sheet, rand(3:length(sheet)-3))
        @test Γ ≈ sum(circulation, (sheet, sheet₊))

        @test Sheets.compute_trapezoidal_weights(sheet.Ss) ≈ circulation.(sheet.blobs)
        @test Sheets.compute_trapezoidal_weights(sheet₊.Ss) ≈ circulation.(sheet₊.blobs)

        N = 200
        θ = range(π, 0, length = N)
        zs = complex.(cos.(θ))
        Γs = sin.(θ)
        Sheets.redistribute_points!(sheet, zs, Γs)
        fill!(Γs, 0.0)
        @test sheet.Ss != Γs
        @test Sheets.compute_trapezoidal_weights(sheet.Ss) ≈ getfield.(sheet.blobs, :S)
        @test Elements.position.(sheet.blobs) == getfield.(sheet.blobs, :z)
        @test sheet.zs == Elements.position.(sheet.blobs)

        Sheets.remesh!(sheet, 0.005)
        @test length(sheet) == 400
        @test sheet.zs ≈ range(-1, 1, length = 400)
        @test norm(sheet.Ss .- sqrt.(1 .- range(-1,1, length=400).^2)) ≤ 1e-3

        θ = range(π, 0, length = N)
        zs = complex.(cos.(θ))
        Γs = sin.(θ)
        Sheets.redistribute_points!(sheet, zs, Γs)
        sheet₁ = deepcopy(sheet)
        sheet₂ = deepcopy(sheet)

        @test_logs (:warn, "Cannot remesh, sheet length smaller than nominal spacing") (:warn, "Filter not applied, total sheet length smaller than nominal spacing") Sheets.filter!(sheet₁, 3.0, 0.03)
        @test sheet₁.zs == sheet₂.zs
        @test sheet₁.Ss == sheet₂.Ss
        @test sheet₁.blobs == sheet₂.blobs

        @test_logs (:warn, "Cannot remesh, sheet length smaller than nominal spacing") Sheets.remesh!(sheet₂, 3.0)
        @test sheet₁.zs == sheet₂.zs
        @test sheet₁.Ss == sheet₂.Ss
        @test sheet₁.blobs == sheet₂.blobs

        Sheets.remesh!(sheet₂, 0.01)
        Sheets.filter_position!(sheet₂, 0.03)

        property = copy(sheet₁.Ss)
        Sheets.filter!(sheet₁, 0.01, 0.03, (property,))

        @test sheet₁.zs == sheet₂.zs
        @test sheet₁.Ss == sheet₂.Ss
        @test sheet₁.Ss == property
    end
end
