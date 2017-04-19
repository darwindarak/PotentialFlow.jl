@testset "Vortex Sheets" begin

    N = 100
    zs = rand(Complex128, N)
    Γs = accumulate(+, rand(N))
    δ = rand()

    sheet = Vortex.Sheet(zs, Γs, δ)

    @test Vortex.circulation(sheet) ≈ Vortex.circulation(sheet.blobs)

end
