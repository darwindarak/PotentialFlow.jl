@testset "Vortex Sheets" begin

    N = 100
    zs = rand(Complex128, N)
    Γs = accumulate(+, rand(N))
    δ = rand()

    sheet = Vortex.Sheet(zs, Γs, δ)

    @test Vortex.circulation(sheet) ≈ Vortex.circulation(sheet.blobs)

    truncate_n = rand(1:50)
    ΣΓ = Vortex.circulation(sheet)
    ΔΓ₀ = Γs[truncate_n] - Γs[1]
    ΔΓ  = Vortex.Sheets.truncate!(sheet, truncate_n)
    @test ΔΓ₀ ≈ ΔΓ

    @test Vortex.Sheets.compute_trapezoidal_weights(sheet.Γs) ≈ getfield.(sheet.blobs, :Γ)
    @test ΔΓ + Vortex.circulation(sheet) ≈ ΣΓ

    ws = allocate_velocity(sheet)
    wb = allocate_velocity(sheet.blobs)

    self_induce_velocity!(ws, sheet)
    for (i, t) in enumerate(sheet.blobs), s in sheet.blobs
        wb[i] += induce_velocity(t, s)
    end
    @test wb ≈ ws
end
