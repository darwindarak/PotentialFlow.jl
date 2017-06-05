@testset "Plate" begin
    @testset "Singular Interactions" begin
        N = 100
        zs = rand(Complex128, N)
        δ  = rand()

        sheet = Vortex.Sheet(zs, cumsum(rand(N)), δ)
        Γs = [b.Γ for b in sheet.blobs]

        points = Vortex.Point.(zs, Γs)
        blobs = Vortex.Blob.(zs, Γs, δ)

        plate = Vortex.Plate(128, 2.0, rand(Complex128), 0.5π*rand())

        vel_s = zeros(Complex128, plate.N)
        vel_p = zeros(Complex128, plate.N)

        Vortex.induce_velocity!(vel_s, plate, sheet)

        Vortex.induce_velocity!(vel_p, plate, points)
        vel_b = Vortex.induce_velocity(plate, blobs)

        @test vel_s == vel_b
        @test vel_p == vel_b

        reset_velocity!(vel_p)
        for (i, b) in enumerate(blobs)
            induce_velocity!(vel_p, plate, b)
        end

        motion = Vortex.Plates.PlateMotion(0.0, 0.0)
        Vortex.Plates.enforce_no_flow_through!(plate, motion, blobs)
        @test vel_p == vel_b
        @test Vortex.Plates.Chebyshev.firstkind(plate.C, plate.ss) ≈ exp(-im*plate.α).*vel_p
    end
    @testset "Bound Circulation" begin
        include("utils/circle_plane.jl")

        c = rand(Complex128)
        ċ = rand(Complex128)
        α = 0.5π*rand()
        α̇ = rand()

        N = 10
        ζs = (2 .+ rand(N)).*exp.(2π.*rand(N))
        zs = c .+ 0.5.*exp(im*α).*(ζs .+ 1./ζs)
        Γs = 1 .- 2.*rand(N)

        Np = 129
        J = JoukowskyMap(c, α)

        Δż_circle = map(linspace(π, 0, Np)) do θ
            η₊ = exp.(im*θ)
            η₋ = conj.(η₊)

            ż₋ = W_vortex(η₋, ζs, Γs) + W_motion(η₋, ċ, α, α̇)
            ż₊ = W_vortex(η₊, ζs, Γs) + W_motion(η₊, ċ, α, α̇)

            ż₊ = conj(ż₊/J(η₊, 1))
            ż₋ = conj(ż₋/J(η₋, 1))

            exp(-im*α).*(ż₋ .- ż₊)
        end;
        @test norm(imag.(Δż_circle[2:end-1])) ≤ 1e-10

        points = Vortex.Point.(zs, Γs)
        plate = Vortex.Plate(Np, 2.0, c, α)
        plate_vel = Vortex.Plates.PlateMotion(ċ, α̇)
        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, points)
        γs = zeros(Float64, Np)
        Vortex.Plates.strength!(γs, plate)
        @test maximum(abs2.(γs[2:end-1] .- real.(Δż_circle[2:end-1]))) ≤ 256eps()

        @test γs == Vortex.Plates.strength(plate)

        kutta_points = Vortex.Point.(c .+ [-1.1, 1.1].*exp(im*α), 1.0)
        Vortex.Plates.vorticity_flux!(plate, kutta_points[1], kutta_points[2])

        ss = linspace(-1, 1, 4001)
        γs = Vortex.Plates.strength(plate, ss)
        Γs = Vortex.Plates.bound_circulation(plate, ss)
        @test Γs[ceil(Int,length(ss)/2)] ≈ Vortex.Plates.bound_circulation(plate)[ceil(Int,length(plate)/2)]

        @test norm(γs - gradient(Γs, step(ss)))/length(ss) < 1e-3
    end

    @testset "Induced Velocities" begin
        c = rand(Complex128)
        ċ = rand(Complex128)
        α = 0.5π*rand()
        α̇ = rand()

        @test c*exp(-im*α) == Vortex.Plates.tangent(c, α) + im*Vortex.Plates.normal(c, α)

        J = JoukowskyMap(c, α)
        N = 100
        ζs = (2 .+ rand(N)).*exp.(2π.*rand(N))
        zs = J.(ζs)
        Γs = 1 .- 2.*rand(N)

        Np = 128
        Nt = 256

        ζt = 10.0.*exp.(im.*linspace(0, 2π, Nt))
        zt = J.(ζt)

        ż_circle = map(ζt) do ζ
            ż = W_vortex(ζ, ζs, Γs) + W_motion(ζ, ċ, α, α̇)
            conj(ż/J(ζ, 1))
        end

        points = Vortex.Point.(zs, Γs)
        plate = Vortex.Plate(Np, 2.0, c, α)
        plate_vel = Vortex.Plates.PlateMotion(ċ, α̇)
        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, points)

        sys = (plate, points)
        żs = Vortex.induce_velocity(zt, sys)
        @test ż_circle ≈ żs
    end

    @testset "Suction Parameters" begin
        U = rand()
        α = rand()*0.5π
        L = 2rand()

        plate = Vortex.Plate(128, L, 0.0, α)
        plate_vel = Vortex.Plates.PlateMotion(U, 0.0)
        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, ())

        C = deepcopy(plate.C)

        point = Vortex.Point(-Inf, 1.0);
        _, Γ, _, _ = Vortex.Plates.vorticity_flux!(plate, point, point, Inf, 0)

        ∂A = Vortex.Plates.influence_on_plate(plate, point)
        C .+= Γ.*∂A

        @test Γ ≈ -π*U*L*sin(α)
        @test Γ == -plate.Γ

        _, b₋ = Vortex.Plates.suction_parameters(plate)
        @test abs(b₋) ≤ eps()
        @test Vortex.Plates.strength(plate, -1) == 0
        @test C ≈ plate.C

        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, ())
        @test_throws AssertionError Vortex.Plates.vorticity_flux(plate, point, point, 0, 0)

        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, ())
        Γ, _, _, _ = Vortex.Plates.vorticity_flux(plate, point, point, 0, Inf)
        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, Vortex.Point(-Inf, Γ))
        b₊, _ = Vortex.Plates.suction_parameters(plate)
        @test abs(b₊) ≤ eps()
        @test Vortex.Plates.strength(plate, 1) == 0

        Vortex.Plates.enforce_no_flow_through!(plate, plate_vel, ())
        Γ₊, Γ₋, _, _ = Vortex.Plates.vorticity_flux(plate, point, point, Inf, Inf)
        @test Γ₊ == Γ₋ == 0
    end

@testset "Impulse" begin
    ċ = rand(Complex128)
    α = rand()*0.5π
    L = 2rand()
    plate = Vortex.Plate(128, L, 0.0, α)
    motion = allocate_velocity(plate)

    motion.ċ = ċ
    motion.α̇ = 0.0

    points = Vortex.Point.(2.0im + rand(Complex128, 20), rand(20))
    sheet  = Vortex.Sheet(2.0im + rand(Complex128, 20), cumsum(rand(20)), rand())

    impulse  = sum(points) do p
        p.Γ*Vortex.unit_impulse(p, plate)
    end

    impulse += sum(sheet.blobs) do b
        b.Γ*Vortex.unit_impulse(b, plate)
    end

    impulse += im*π*0.5L*imag(ċ*exp(-im*α))
    impulse *= 0.5L*exp(im*α)

    Vortex.Plates.enforce_no_flow_through!(plate, motion, (points, sheet))

    @test Vortex.circulation((plate, points, sheet)) ≤ 128eps()
    @test norm(impulse .- Vortex.impulse((plate, points, sheet))) ≤ 1e-5
end

@testset "Advection" begin
    motion = Vortex.Plates.PlateMotion(rand(Complex128), rand())
    c = rand(Complex128)
    α = 0.5π*rand()
    L = 2rand()

    points = Vortex.Point.(rand(Complex128, 10), rand(10))

    plate  = Vortex.Plate(128, L, c, α)
    Vortex.Plates.enforce_no_flow_through!(plate, motion, points)

    # Internal variables should be overwritten/resized on `advect!`
    plate₊ = Vortex.Plate(10, rand(), rand(), rand())

    Δt = 1e-2
    advect!(plate₊, plate, motion, Δt)
    advect!(plate,  plate, motion, Δt)

    @test plate.zs ≈ @. plate.c + 0.5plate.L*exp(im*plate.α)*plate.ss

    # Check that interval values are the same
    @test length(plate₊) == length(plate)
    @test plate₊.C  == plate.C
    @test plate₊.zs == plate.zs
    @test plate₊.c  == plate.c
    @test plate₊.α  == plate.α
    @test plate₊.Γ  == plate.Γ
    @test plate₊.B₀ == plate.B₀
    @test plate₊.B₁ == plate.B₁

    # Check that internal pointers are different
    @test plate₊.C  ≢ plate.C
    @test plate₊.zs ≢ plate.zs
    @test plate₊.ss ≢ plate.ss
    @test plate₊.A.data ≡ plate₊.C
end

@testset "[De]serialization" begin
    motion = Vortex.Plates.PlateMotion(rand(Complex128), rand())
    c = rand(Complex128)
    α = 0.5π*rand()
    L = 2rand()
    plate  = Vortex.Plate(128, L, c, α)
    Vortex.Plates.enforce_no_flow_through!(plate, motion, ())

    buff = IOBuffer()
    serialize(buff, plate)
    seekstart(buff)

    plate₊ = deserialize(buff)
    Vortex.Plates.enforce_no_flow_through!(plate₊, motion, ())
    @test plate₊.C ≈ plate.C
end

end
