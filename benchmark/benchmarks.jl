using VortexModel

@benchgroup "Pairwise-dispatch" begin
    N = 500
    points = Vortex.Point.(rand(Complex128, N), rand(N))
    sheet  = Vortex.Sheet(rand(Complex128, N), accumulate(+, rand(N)), rand())
    blobs = sheet.blobs

    z = rand(Complex128)

    @bench "points to position" Vortex.induce_velocity(z, points)
    @bench "group to position" Vortex.induce_velocity(z, (points, blobs, sheet))
    @bench "nested group to position" Vortex.induce_velocity(z, (points, blobs, sheet))

    zs = getfield.(points, :z)
    @bench "points to positions" Vortex.induce_velocity(zs, points)
    @bench "sheet to positions" Vortex.induce_velocity(zs, sheet)

    wss = Vortex.allocate_velocity(points)
    @bench "in-place sheet to points" Vortex.induce_velocity!(wss, points, sheet)

    @bench "in-place group to points" Vortex.induce_velocity!(wss, points, (blobs, sheet))

    sys = (sheet, blobs)
    ws = Vortex.allocate_velocity(sys)
    @bench "in-place points to group" Vortex.induce_velocity!(ws, sys, points)

    @bench "self-induced velocities" Vortex.self_induce_velocity!(ws, sys)
end
