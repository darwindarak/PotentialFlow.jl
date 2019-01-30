@testset "Body" begin
    @testset "Geometry" begin
    a1 = 1; b1 = 0.1
    ccoeff = ComplexF64[0.5(a1+b1),0,0.5(a1-b1)]
    b = Bodies.ConformalBody(ccoeff)
    ζ = exp(im*0)
    v = exp(im*π/4)
    @test Bodies.normal(ζ,v,b) ≈ 0.7071067811865475
    @test Bodies.tangent(ζ,v,b) ≈ 0.7071067811865475
    @test conftransform(ζ,b) ≈ 1.0
    p = Bodies.Polygon([-1.0,1.0,1.0,-1.0],[-1.0,-1.0,1.0,1.0])
    b = Bodies.ConformalBody(p)
    ζ = exp(im*0)
    v = exp(im*π/4)
    @test Bodies.normal(ζ,v,b) ≈ 0.7071067811865475
    @test Bodies.tangent(ζ,v,b) ≈ 0.7071067811865475
    @test conftransform(ζ,b) ≈ 1.0
    end

    @testset "Motion" begin
    a1 = 1; b1 = 0.1
    ccoeff = ComplexF64[0.5(a1+b1),0,0.5(a1-b1)]
    b = Bodies.ConformalBody(ccoeff)
    Δt = 1.0
    ċ = 1.0im
    motion = RigidBodyMotion(ċ, 0.0);
    b₊ = deepcopy(b)
    advect!(b₊,b,motion,Δt);
    @test b₊.c ≈ 0.0 + 1.0im
    @test b₊.α ≈ 0.0
    ċ = 0.0im
    motion = RigidBodyMotion(ċ, 1.0);
    b₊ = deepcopy(b)
    advect!(b₊,b,motion,Δt);
    @test b₊.c ≈ 0.0 + 0.0im
    @test b₊.α ≈ 1.0
    end


end
