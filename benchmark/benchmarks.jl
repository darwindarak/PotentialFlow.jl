using PotentialFlow

@benchgroup "Pairwise-dispatch" begin
    N = 500
    t = 0.0

    points = Vortex.Point.(rand(ComplexF64, N), rand(N))
    sheet  = Vortex.Sheet(rand(ComplexF64, N), accumulate(+, rand(N)), rand())
    blobs = sheet.blobs

    z = rand(ComplexF64)

    @bench "points to position" induce_velocity(z, points, t)
    @bench "group to position" induce_velocity(z, (points, blobs, sheet), t)
    @bench "nested group to position" induce_velocity(z, (points, blobs, sheet), t)

    zs = getfield.(points, :z)
    @bench "points to positions" induce_velocity(zs, points, t)
    @bench "sheet to positions" induce_velocity(zs, sheet, t)

    wss = allocate_velocity(points)
    @bench "in-place sheet to points" induce_velocity!(wss, points, sheet, t)

    @bench "in-place group to points" induce_velocity!(wss, points, (blobs, sheet), t)

    sys = (sheet, blobs)
    ws = allocate_velocity(sys)
    @bench "in-place points to group" induce_velocity!(ws, sys, points, t)

    @bench "self-induced velocities" self_induce_velocity!(ws, sys, t)
end
