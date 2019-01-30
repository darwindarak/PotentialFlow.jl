@testset "Pairwise Dispatch" begin
    N = 100
    points = Vortex.Point.(rand(ComplexF64, N), rand(N))
    sheet  = Vortex.Sheet(rand(ComplexF64, N), accumulate(+, rand(N)), rand())
    blobs = sheet.blobs

    z = rand(ComplexF64)

    wp = induce_velocity(z, points, 0.0)
    wb = induce_velocity(z, blobs, 0.0)
    ws = induce_velocity(z, sheet, 0.0)

    @test wb == ws

    w = induce_velocity(z, (points, blobs, sheet), 0.0)
    @test w == (wp + wb + ws)

    w = induce_velocity(z, ((points, blobs), sheet), 0.0)
    @test w == (wp + wb + ws)

    zs = getfield.(points, :z)
    wps = induce_velocity(zs, points, 0.0)
    wbs = induce_velocity(zs, blobs, 0.0)
    wss = induce_velocity(zs, sheet, 0.0)

    @test wss == wbs
    @test wbs == induce_velocity(points, blobs, 0.0)
    @test wss == induce_velocity(points, sheet, 0.0)

    @test wbs .+ wss == induce_velocity(points, (blobs, sheet), 0.0)

    reset_velocity!(wss)
    induce_velocity!(wss, points, sheet, 0.0)

    @test wss == wbs

    reset_velocity!(wps)
    induce_velocity!(wps, points, (blobs, sheet), 0.0)
    @test wps == wss .+ wbs

    # Distribute to a group of sources
    sys = (sheet, blobs)
    ws = allocate_velocity(sys)
    induce_velocity!(ws, sys, points, 0.0)

    @test ws[1] == induce_velocity(sheet, points, 0.0)
    @test ws[2] == induce_velocity(blobs, points, 0.0)

    # Self-induce velocities
    reset_velocity!(ws)
    self_induce_velocity!(ws, sys, 0.0)
    @test ws == self_induce_velocity(sys, 0.0)

    wb = allocate_velocity(blobs)
    induce_velocity!(wb, blobs, blobs, 0.0)
    induce_velocity!(wb, blobs, sheet, 0.0)

    @test ws[1] ≈ wb

    reset_velocity!(wb)
    induce_velocity!(wb, sheet, sheet, 0.0)
    induce_velocity!(wb, sheet, blobs, 0.0)

    @test ws[2] ≈ wb

    points = Vortex.Point.(rand(ComplexF64, 2N), rand(2N))
    reset_velocity!(ws, (points, blobs))
    @test length.(ws) == (2N, N)

    reset_velocity!(ws, (points, points))
    reset_velocity!(wb, points)
    self_induce_velocity!(wb, points, 0.0)
    mutually_induce_velocity!(ws[1], ws[2], points, points, 0.0)
    @test ws[1] ≈ ws[2]
    @test ws[1] ≈ wb

    reset_velocity!(wb, points)
    induce_velocity!(wb, points, points, 0.0)
    @test ws[1] ≈ wb

    @testset "defaults" begin
        # Minimum implementation of a vortex point source
        @eval struct NewPoint <: Elements.Element
            z::ComplexF64
            Γ::Float64
        end

        Elements.position(p::NewPoint) = p.z
        Elements.circulation(p::NewPoint) = p.Γ
        Elements.impulse(p::NewPoint) = -im*p.z*p.Γ

        Elements.kind(::NewPoint) = Elements.Singleton
        Elements.kind(::Type{NewPoint}) = Elements.Singleton

        function Motions.induce_velocity(z::ComplexF64, p::NewPoint, t)
            p.Γ*Vortex.Points.cauchy_kernel(z - p.z)
        end

        zs = rand(ComplexF64, 100)
        Γs = rand(100)

        new_points = NewPoint.(zs, Γs)
        points = Vortex.Point.(zs, Γs)

        wn = allocate_velocity(new_points)
        wp = allocate_velocity(points)

        self_induce_velocity!(wn, new_points, 0.0)
        self_induce_velocity!(wp, points, 0.0)

        @test wn ≈ wp
    end

    @testset "Property Macro" begin
        import PotentialFlow.Properties: @property

        @eval Elements (@property begin
                      signature = max_blob_radius(s::Source)
                      stype = Int
                      reduce = max
                      end)

        @eval Elements (@property begin
                      signature = induce_count(t::Target, tc::Target, s::Source, sc::Source)
                      preallocator = allocate_count
                      stype = Int
                      end)

        @eval Vortex.Points begin
            Elements.induce_count(::ComplexF64, tcount, ::Point, scount) = scount
            Elements.max_blob_radius(::Point) = 0
        end
        @eval Vortex.Blobs begin
            Elements.induce_count(::ComplexF64, tcount, ::Blob, scount) = tcount
            Elements.max_blob_radius(b::Blob) = b.δ
        end

        N = rand(1:100)
        points = Vortex.Point.(rand(ComplexF64, N), rand(N))
        pcounts = fill(1, length(points))

        blobs  = Vortex.Blob.(rand(ComplexF64, N), rand(N), rand())
        bcounts = fill(2, length(blobs))

        sources  = Source.Blob.(rand(ComplexF64, N), rand(N), rand())
        bcounts = fill(2, length(blobs))

        targets = rand(ComplexF64, rand(1:100))
        tcounts = fill(3, length(targets))

        @test Elements.max_blob_radius(points) == 0
        @test Elements.max_blob_radius(blobs) == blobs[1].δ
        @test Elements.max_blob_radius((points, blobs)) == blobs[1].δ

        @test Elements.induce_count(targets[1], tcounts[1], points, pcounts) == length(points)
        @test Elements.induce_count(targets[1], tcounts[1], blobs, bcounts) == 3length(blobs)
        @test Elements.induce_count(targets[1], tcounts[1],
                                  (points, blobs), (pcounts, bcounts)) == 3length(blobs) + length(points)

        out = Elements.induce_count(targets, tcounts, (points, blobs), (pcounts, bcounts))
        @test length(out) == length(targets)

        out₂ = Elements.allocate_count(targets)
        Elements.induce_count!(out₂, targets, tcounts, (points, blobs), (pcounts, bcounts))

        @test out == out₂
        @test all(out .== 3length(blobs) + length(points))
    end
end
