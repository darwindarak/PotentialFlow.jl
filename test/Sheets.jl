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
end
