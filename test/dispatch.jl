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
    @test ws == Vortex.self_induce_velocity(sys)

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

    Vortex.reset_velocity!(ws, (points, points))
    Vortex.reset_velocity!(wb, points)
    Vortex.self_induce_velocity!(wb, points)
    Vortex.mutually_induce_velocity!(ws[1], ws[2], points, points)
    @test ws[1] ≈ ws[2]
    @test ws[1] ≈ wb

    Vortex.reset_velocity!(wb, points)
    Vortex.induce_velocity!(wb, points, points)
    @test ws[1] ≈ wb

    @testset "defaults" begin
        # Minimum implementation of a vortex point source
        @eval struct NewPoint <: Vortex.Element
            z::Complex128
            Γ::Float64
        end

        Vortex.position(p::NewPoint) = p.z
        Vortex.circulation(p::NewPoint) = p.Γ
        Vortex.impulse(p::NewPoint) = -im*p.z*p.Γ

        Vortex.kind(::NewPoint) = Vortex.Singleton
        Vortex.kind(::Type{NewPoint}) = Vortex.Singleton

        function Vortex.induce_velocity(z::Complex128, p::NewPoint)
            p.Γ*Vortex.Points.cauchy_kernel(z - p.z)
        end

        zs = rand(Complex128, 100)
        Γs = rand(100)

        new_points = NewPoint.(zs, Γs)
        points = Vortex.Point.(zs, Γs)

        wn = allocate_velocity(new_points)
        wp = allocate_velocity(points)

        self_induce_velocity!(wn, new_points)
        self_induce_velocity!(wp, points)

        @test wn ≈ wp
    end

    @testset "Property Macro" begin
        @test_throws ArgumentError @eval Vortex (@property prop count Int)

        @eval Vortex (@property induced count Int)
        @eval Vortex.Points begin
            Vortex.induce_count(::Complex128, ::Point) = 1
        end
        @eval Vortex.Blobs begin
            Vortex.induce_count(::Complex128, ::Blob) = 1
        end

        N = rand(1:100)
        points = Vortex.Point.(rand(Complex128, N), rand(N))
        blobs  = Vortex.Blob.(rand(Complex128, N), rand(N), rand())

        targets = rand(Complex128, rand(1:100))

        @test Vortex.induce_count(targets[1], points) == N
        @test Vortex.induce_count(targets[1], blobs) == N
        @test Vortex.induce_count(targets[1], (points, blobs)) == 2N

        out = Vortex.induce_count(targets, (points, blobs))
        @test length(out) == length(targets)

        out₂ = Vortex.allocate_count(targets)
        Vortex.induce_count!(out₂, targets, (points, blobs))

        @test out == out₂
        @test all(out .== 2N)
    end
end
