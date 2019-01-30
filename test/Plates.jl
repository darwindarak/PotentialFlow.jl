@testset "Plate" begin

    @testset "Singular Interactions" begin
        N = 100
        zs = rand(ComplexF64, N)
        δ  = rand()

        sheet = Vortex.Sheet(zs, cumsum(rand(N)), δ)
        Γs = circulation.(sheet.blobs)

        points = Vortex.Point.(zs, Γs)
        blobs = Vortex.Blob.(zs, Γs, δ)

        plate = Plate(128, 2.0, rand(ComplexF64), 0.5π*rand())

        vel_s = zeros(ComplexF64, plate.N)
        vel_p = zeros(ComplexF64, plate.N)

        induce_velocity!(vel_s, plate, sheet, 0)

        induce_velocity!(vel_p, plate, points, 0)
        vel_b = induce_velocity(plate, blobs, 0)

        @test vel_s == vel_b
        @test vel_p == vel_b

        reset_velocity!(vel_p)
        for (i, b) in enumerate(blobs)
            induce_velocity!(vel_p, plate, b, 0)
        end

        motion = Plates.RigidBodyMotion(0.0, 0.0)
        Plates.enforce_no_flow_through!(plate, motion, blobs, 0)
        @test vel_p == vel_b
        @test Plates.Chebyshev.firstkind(plate.C, plate.ss) ≈ exp(-im*plate.α).*vel_p
    end
    @testset "Bound Circulation" begin
      import PotentialFlow.Utils:centraldiff

        include("utils/circle_plane.jl")


        c = rand(ComplexF64)
        ċ = rand(ComplexF64)
        α = 0.5π*rand()
        α̇ = rand()

        N = 10
        ζs = (2 .+ rand(N)).*exp.(2π.*rand(N))
        zs = c .+ 0.5.*exp(im*α).*(ζs .+ 1 ./ ζs)
        Γs = 1 .- 2 .* rand(N)

        Np = 129
        J = JoukowskyMap(c, α)

        Δż_circle = map(range(π, 0, length = Np)) do θ
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
        plate = Plate(Np, 2.0, c, α)
        plate_vel = Plates.RigidBodyMotion(ċ, α̇)
        Plates.enforce_no_flow_through!(plate, plate_vel, points, 0)
        γs = zeros(Float64, Np)
        Plates.strength!(γs, plate)
        @test maximum(abs2.(γs[2:end-1] .- real.(Δż_circle[2:end-1]))) ≤ 256eps()

        @test γs == Plates.strength(plate)

        kutta_points = Vortex.Point.(c .+ [-1.1, 1.1].*exp(im*α), 1.0)
        Plates.vorticity_flux!(plate, kutta_points[1], kutta_points[2], 0)

        ss = range(-1, 1, length = 4001)
        γs = Plates.strength(plate, ss)
        Γs = Plates.bound_circulation(plate, ss)
        @test Γs[ceil(Int,length(ss)/2)] ≈ Plates.bound_circulation(plate)[ceil(Int,length(plate)/2)]

        @test norm(γs - centraldiff(Γs,step(ss)))/length(ss) < 1e-3
    end

    @testset "Induced Velocities" begin
        c = rand(ComplexF64)
        ċ = rand(ComplexF64)
        α = 0.5π*rand()
        α̇ = rand()

        @test c*exp(-im*α) == Plates.tangent(c, α) + im*Plates.normal(c, α)

        J = JoukowskyMap(c, α)
        N = 100
        ζs = (2 .+ rand(N)).*exp.(2π.*rand(N))
        zs = J.(ζs)
        Γs = 1 .- 2 .* rand(N)

        Np = 128
        Nt = 256

        ζt = 10.0.*exp.(im.*range(0, 2π, length = Nt))
        zt = J.(ζt)

        ż_circle = map(ζt) do ζ
            ż = W_vortex(ζ, ζs, Γs) + W_motion(ζ, ċ, α, α̇)
            conj(ż/J(ζ, 1))
        end

        points = Vortex.Point.(zs, Γs)
        plate = Plate(Np, 2.0, c, α)
        plate_vel = Plates.RigidBodyMotion(ċ, α̇)
        Plates.enforce_no_flow_through!(plate, plate_vel, points, 0)

        sys = (plate, points)
        żs = induce_velocity(zt, sys, 0)
        @test ż_circle ≈ żs
    end

    @testset "Suction Parameters" begin
        U = rand()
        α = rand()*0.5π
        L = 2rand()

        plate = Plate(128, L, 0.0, α)
        plate_vel = Plates.RigidBodyMotion(U, 0.0)
        Plates.enforce_no_flow_through!(plate, plate_vel, (), 0)

        C = deepcopy(plate.C)

        point = Vortex.Point(-Inf, 1.0);
        _, Γ, _, _ = Plates.vorticity_flux!(plate, point, point, 0, Inf, 0)

        ∂A = Plates.influence_on_plate(plate, point, 0)
        C .+= Γ.*∂A

        @test Γ ≈ -π*U*L*sin(α)
        @test Γ == -plate.Γ

        _, b₋ = Plates.suction_parameters(plate)
        @test abs(b₋) ≤ eps()
        @test Plates.strength(plate, -1) == 0
        @test C ≈ plate.C

        Plates.enforce_no_flow_through!(plate, plate_vel, (), 0)
        @test_throws AssertionError Plates.vorticity_flux(plate, point, point, 0, 0, 0)

        Plates.enforce_no_flow_through!(plate, plate_vel, (), 0)
        Γ, _, _, _ = Plates.vorticity_flux(plate, point, point, 0, 0, Inf)
        Plates.enforce_no_flow_through!(plate, plate_vel, Vortex.Point(-Inf, Γ), 0)
        b₊, _ = Plates.suction_parameters(plate)
        @test abs(b₊) ≤ eps()
        @test Plates.strength(plate, 1) == 0

        Plates.enforce_no_flow_through!(plate, plate_vel, (), 0)
        Γ₊, Γ₋, _, _ = Plates.vorticity_flux(plate, point, point, 0, Inf, Inf)
        @test Γ₊ == Γ₋ == 0
    end

@testset "Impulse" begin
    ċ = rand(ComplexF64)
    α = rand()*0.5π
    L = 2rand()
    plate = Plate(128, L, 0.0, α)
    @test_logs (:warn,"Plate kinematics should be initialized manually.  This simply returns a stationary motion") allocate_velocity(plate)

    motion = Plates.RigidBodyMotion(ċ, 0.0)

    points = Vortex.Point.(2.0im .+ rand(ComplexF64, 20), rand(20))
    sheet  = Vortex.Sheet(2.0im .+ rand(ComplexF64, 20), cumsum(rand(20)), rand())

    impulse  = sum(points) do p
        circulation(p)*Plates.unit_impulse(p, plate)
    end

    impulse += sum(sheet.blobs) do b
        circulation(b)*Plates.unit_impulse(b, plate)
    end

    impulse += im*π*0.5L*imag(ċ*exp(-im*α))
    impulse *= 0.5L*exp(im*α)

    Plates.enforce_no_flow_through!(plate, motion, (points, sheet), 0)

    @test circulation((plate, points, sheet)) ≤ 128eps()
    @test norm(impulse .- Plates.impulse((plate, points, sheet))) ≤ 1e-5
end

@testset "Advection" begin
    motion = Plates.RigidBodyMotion(rand(ComplexF64), rand())
    c = rand(ComplexF64)
    α = 0.5π*rand()
    L = 2rand()

    points = Vortex.Point.(rand(ComplexF64, 10), rand(10))

    plate = Plate(128, L, c, α)
    Plates.enforce_no_flow_through!(plate, motion, points, 0)

    # Internal variables should be overwritten/resized on `advect!`
    plate₊ = Plate(10, rand(), rand(), rand())

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

end
