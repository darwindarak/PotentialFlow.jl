@testset "Pairwise Dispatch" begin
    N = 100
    points = Vortex.Point.(rand(Complex128, N), rand(N))
    sheet  = Vortex.Sheet(rand(Complex128, N), accumulate(+, rand(N)), rand())
    blobs = sheet.blobs

    z = rand(Complex128)

    wp = Vortex.induce_velocity(z, points)
    wb = Vortex.induce_velocity(z, blobs)
    ws = Vortex.induce_velocity(z, sheet)

    @test wb == ws

    w = Vortex.induce_velocity(z, (points, blobs, sheet))
    @test w == (wp + wb + ws)

    w = Vortex.induce_velocity(z, ((points, blobs), sheet))
    @test w == (wp + wb + ws)

    zs = getfield.(points, :z)
    wps = Vortex.induce_velocity(zs, points)
    wbs = Vortex.induce_velocity(zs, blobs)
    wss = Vortex.induce_velocity(zs, sheet)

    @test wss == wbs
    @test wbs == Vortex.induce_velocity(points, blobs)
    @test wss == Vortex.induce_velocity(points, sheet)

    @test wbs .+ wss == Vortex.induce_velocity(points, (blobs, sheet))

    Vortex.reset_velocity!(wss)
    Vortex.induce_velocity!(wss, points, sheet)

    @test wss == wbs

    Vortex.reset_velocity!(wps)
    Vortex.induce_velocity!(wps, points, (blobs, sheet))
    @test wps == wss .+ wbs
    
    # Distribute to a group of sources
    sys = (sheet, blobs)
    ws = Vortex.allocate_velocity(sys)
    Vortex.induce_velocity!(ws, sys, points)

    @test ws[1] == Vortex.induce_velocity(sheet, points)
    @test ws[2] == Vortex.induce_velocity(blobs, points)

    # Self-induce velocities
    Vortex.reset_velocity!(ws)
    Vortex.self_induce_velocity!(ws, sys)

    wb = Vortex.allocate_velocity(blobs)
    Vortex.induce_velocity!(wb, blobs, blobs)
    Vortex.induce_velocity!(wb, blobs, sheet)

    @test ws[1] ≈ wb

    Vortex.reset_velocity!(wb)
    Vortex.induce_velocity!(wb, sheet, sheet)
    Vortex.induce_velocity!(wb, sheet, blobs)

    @test ws[2] ≈ wb

    points = Vortex.Point.(rand(Complex128, 2N), rand(2N))
    Vortex.reset_velocity!(ws, (points, blobs))
    @test length.(ws) == (2N, N)

end
