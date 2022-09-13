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

    @testset "Derivatives" begin

    a1 = 0.5; b1 = 0.1; ccoeff = ComplexF64[0.5(a1+b1),0,0.5(a1-b1)]
    b = Bodies.ConformalBody(ccoeff,ComplexF64(0.0),π/4)

    vort_z = [Vortex.Blob(0.5 + 0.2im, 1.0, 0.01),
              Vortex.Blob(0.2im, 1.0, 0.01),
              Vortex.Blob(-1.0im, -2.0, 0.01)]

    vort_ζ = Elements.inverse_conftransform(vort_z,b)

    ζ = 1.0+0.1im
    @test Bodies.dpdzv(ζ,1,vort_ζ,b) ≈ -1.3941088515895026 + 1.1089339580694477im
    @test Bodies.dpdzv(ζ,2,vort_ζ,b) ≈ 0.8402098734913479 + 0.8115186077355655im

    motion = RigidBodyMotion(RigidBodyMotions.UnsteadyTransRot(-0.5 + 1.0im,0.0im,0.2,0.0))
    Bodies.enforce_no_flow_through!(b, motion, Element[], 0.0)
    #@test Bodies.dpdzv(ζ,1,vort_ζ,b) ≈ -14.91887479595418 + 11.326836023135279im
    #@test Bodies.dpdΓv(ζ,1,vort_ζ,b) ≈ -0.10685330645316929
    @test Bodies.dpdzv(ζ,1,vort_ζ,b) ≈ -13.564956302107685 + 18.2006953682202im
    @test Bodies.dpdΓv(ζ,1,vort_ζ,b) ≈ -10.74045898676263

    #@test Bodies.dpdU(ζ,1,vort_ζ,b) ≈ -7.445770018204262
    #@test Bodies.dpdU(ζ,2,vort_ζ,b) ≈ 7.489634497871336
    #@test Bodies.dpdU(ζ,3,vort_ζ,b) ≈ -22.042738470789313
    @test Bodies.dpdU(ζ,1,vort_ζ,b) ≈ -8.548018373889507
    @test Bodies.dpdU(ζ,2,vort_ζ,b) ≈ 4.126900666032036
    @test Bodies.dpdU(ζ,3,vort_ζ,b) ≈ -30.61646732613466

    @test Bodies.dpdUdot(ζ,1,vort_ζ,b) ≈ 0.01176355259288305
    @test Bodies.dpdUdot(ζ,2,vort_ζ,b) ≈ 0.09900990099009899
    @test Bodies.dpdUdot(ζ,3,vort_ζ,b) ≈ 0.04950495049504951


    end

    @testset "Force and moment" begin
    a1 = 0.5; b1 = 0.1; ccoeff = ComplexF64[0.5(a1+b1),0,0.5(a1-b1)]
    b = Bodies.ConformalBody(ccoeff,ComplexF64(0.0),π/4)

    vort_z = [Vortex.Point(0.5 + 0.2im, 1.0),
              Vortex.Point(0.2im, 1.0),
              Vortex.Point(-1.0im, -2.0)]

    vort_ζ = Elements.inverse_conftransform(vort_z,b)

    fx, fy, mr = Bodies.force(vort_ζ,b)

    θ = range(0,2π,length=2001)
    dΘ = θ[2]-θ[1]
    ζc = exp.(im*θ[1:end-1])

    dz̃c = map(ζ -> b.dm(ζ)[1],ζc)
    f_test = -sum(Bodies.pressure(ζc,vort_ζ,b).*dz̃c.*ζc*dΘ)

    z̃c = map(ζ -> b.m(ζ),ζc)
    mr_test = -sum(Bodies.pressure(ζc,vort_ζ,b).*imag(conj(z̃c).*dz̃c.*ζc)*dΘ)

    @test isapprox(fx,real(f_test),atol=1e-8)
    @test isapprox(fy,imag(f_test),atol=1e-8)
    @test isapprox(mr,mr_test,atol=1e-8)

    motion = RigidBodyMotion(RigidBodyMotions.UnsteadyTransRot(-0.2 - 0.3im,0.3im,0.5,-0.2))
    Bodies.enforce_no_flow_through!(b, motion, Element[], 0.0)

    fx, fy, mr = Bodies.force(vort_ζ,b)
    f_test = -sum(Bodies.pressure(ζc,vort_ζ,b).*dz̃c.*ζc*dΘ)
    mr_test = -sum(Bodies.pressure(ζc,vort_ζ,b).*imag(conj(z̃c).*dz̃c.*ζc)*dΘ)

    @test isapprox(fx,real(f_test),atol=1e-8)
    @test isapprox(fy,imag(f_test),atol=1e-8)
    @test isapprox(mr,mr_test,atol=1e-8)

    end

end
