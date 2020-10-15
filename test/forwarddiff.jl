using LinearAlgebra

import PotentialFlow.Utils: derivative, extract_derivative, value, partials,
          Dual,ComplexComplexDual,ComplexRealDual

const DELTA=1e-9
const BIGEPS = 1000*eps(1.0)
const TOL=5e-6
const BIGTOL=1e-5
const BIGGESTTOL=1e-3

safenorm(a) = norm(filter(x -> ~isnan(x),a))

@testset "Complex Automatic Differentiation" begin

  @testset "Basic derivatives" begin

    z = rand(ComplexF64)

    dz1, dzstar1 = derivative(z -> log(sqrt(z)),z)

    dzex = 1/(2*z)

    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1),0,atol=BIGEPS)

    dz1, dzstar1 = derivative(z -> conj(z),z)
    @test dzstar1 == one(z) && dz1 == zero(z)

    f(z) = 0.5im/(π*conj(z))

    dz1, dzstar1 = derivative(f,z)
    dzstarex = -0.5im/(pi*conj(z)^2)

    @test isapprox(abs(dz1),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    z0 = rand(ComplexF64)
    dz1, dzstar1 = derivative(z -> f(z-z0),z)

    dzstarex = -0.5im/(pi*conj(z-z0)^2)

    @test isapprox(abs(dz1),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    dz1, dzstar1 = derivative(abs,z)

    dzex = 0.5*conj(z)/abs(z)
    dzstarex = 0.5*z/abs(z)
    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    dz1, dzstar1 = derivative(z -> z + sqrt(log(z)) - 1/z^2,z)

    dzex = 1 + 1/(2*sqrt(log(z))*z) + 2/z^3
    dzstarex = 0.0

    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    dz1, dzstar1 = derivative(z -> real(z - √(z - 1)*√(z + 1)),z)

    dzex = 0.5*(1 - 0.5*sqrt((z+1)/(z-1)) - 0.5*sqrt((z-1)/(z+1)))
    dzstarex = 0.5*conj(1 - 0.5*sqrt((z+1)/(z-1)) - 0.5*sqrt((z-1)/(z+1)))

    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

  end

  nblob = 5
  pos = rand(ComplexF64,nblob)
  str = rand(length(pos))

  σ = 1e-2
  blobs = Vortex.Blob.(pos,str,σ)
  z = rand(ComplexF64)

  i = rand(1:nblob)

  # For finite difference approximations
  dz = zeros(ComplexF64,length(blobs))
  dz[i] = DELTA
  dΓ = zeros(Float64,length(blobs))
  dΓ[i] = DELTA
  blobsx⁺ = Vortex.Blob.(Elements.position(blobs).+dz,Elements.circulation.(blobs),σ)
  blobsy⁺ = Vortex.Blob.(Elements.position(blobs).+im*dz,Elements.circulation.(blobs),σ)
  blobsΓ⁺ = Vortex.Blob.(Elements.position(blobs),Elements.circulation.(blobs).+dΓ,σ)



  @testset "Basic operations with duals" begin


    dualpos = one(ComplexComplexDual{Nothing},Elements.position(blobs)[i])

    @test value(dualpos) == Elements.position(blobs)[i]
    pr, pi = partials(dualpos)
    @test pr == [1.0,0.0] && pi == [0.0,1.0]

    newblob = Vortex.Blob(dualpos,Elements.circulation(blobs[i]),σ)
    @test value(Elements.position(newblob)) == Elements.position(blobs)[i]

    dwdz,dwdzstar  = extract_derivative(Nothing,induce_velocity(z,newblob,0.0))

    w2 = Elements.circulation(newblob)*PotentialFlow.Blobs.blob_kernel(z - Elements.position(newblob),σ)
    dwdz2, dwdzstar2 = extract_derivative(Nothing,w2)

    @test isapprox(abs(dwdz-dwdz2),0.0,atol=BIGEPS)
    @test isapprox(abs(dwdzstar-dwdzstar2),0.0,atol=BIGEPS)

end

@testset "Induced velocity at a point" begin

    # Finite difference approximations
    dwdx_fd = (induce_velocity(z,blobsx⁺,0.0) - induce_velocity(z,blobs,0.0))/dz[i]
    dwdy_fd =(induce_velocity(z,blobsy⁺,0.0) - induce_velocity(z,blobs,0.0))/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)
    dwdΓ_fd = (induce_velocity(z,blobsΓ⁺,0.0) - induce_velocity(z,blobs,0.0))/dΓ[i]

    # Auto differentation
    newblobs = Vortex.dualize_position(blobs,i,Nothing)
    dwdz, dwdzstar = extract_derivative(Nothing,induce_velocity(z,newblobs,0.0))

    @test isapprox(abs(dwdz-dwdz_fd),0.0,atol=TOL)
    @test isapprox(abs(dwdzstar-dwdzstar_fd),0.0,atol=TOL)

    newblobs = Vortex.dualize_strength(blobs,i,Nothing)

    @test sum(value.(Vortex.circulation.(newblobs))) -
              value(Vortex.circulation(newblobs)) == 0

    dwdΓ = extract_derivative(Nothing,induce_velocity(z,newblobs,0.0))
    @test isapprox(abs(dwdΓ-dwdΓ_fd),0.0,atol=TOL)

end

N = 7
L = 2.0
c = complex(0)
α = 0.0
ċ = complex(0.0)
α̇ = 0.0
p = PotentialFlow.Plate(N,L,c,α)
motion = PotentialFlow.RigidBodyMotion(ċ,α̇)

Plates.enforce_no_flow_through!(p, motion, blobs, 0.0)
px⁺ = deepcopy(p)
Plates.enforce_no_flow_through!(px⁺, motion, blobsx⁺, 0.0)
py⁺ = deepcopy(p)
Plates.enforce_no_flow_through!(py⁺, motion, blobsy⁺, 0.0)
pΓ⁺ = deepcopy(p)
Plates.enforce_no_flow_through!(pΓ⁺, motion, blobsΓ⁺, 0.0)

@testset "Induced velocity at plate points" begin

    C  = zeros(ComplexF64, N)
    dchebt! = Plates.Chebyshev.plan_transform!(C)

    C  = zeros(ComplexF64, N)
    induce_velocity!(C,p,blobs,0.0)

    Cx⁺  = zeros(ComplexF64, N)
    induce_velocity!(Cx⁺,p,blobsx⁺,0.0)
    dwdx_fd = (Cx⁺ - C)/dz[i]

    Cy⁺  = zeros(ComplexF64, N)
    induce_velocity!(Cy⁺,p,blobsy⁺,0.0)
    dwdy_fd = (Cy⁺ - C)/dz[i]

    CΓ⁺  = zeros(ComplexF64, N)
    induce_velocity!(CΓ⁺,p,blobsΓ⁺,0.0)

    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)

    newblobs = Vortex.dualize_position(blobs,i,Nothing)
    C2 = zeros(typeof(ComplexComplexDual()),N)
    induce_velocity!(C2,p,newblobs,0.0)
    dwdz, dwdzstar = extract_derivative(Nothing,C2)

    # test that the induced velocities and their derivatives match
    @test isapprox(norm(value.(C2) - C),0.0,atol=BIGEPS)
    @test isapprox(norm(dwdz - dwdz_fd),0.0,atol=TOL)
    @test isapprox(norm(dwdzstar - dwdzstar_fd),0.0,atol=TOL)

    # finite diff
    dchebt! * C
    dchebt! * Cx⁺
    dchebt! * Cy⁺
    dchebt! * CΓ⁺
    dCdx_fd = (Cx⁺ - C)/dz[i]
    dCdy_fd = (Cy⁺ - C)/dz[i]
    dCdΓ_fd = (CΓ⁺ - C)/dΓ[i]
    dCdz_fd = 0.5*(dCdx_fd - im*dCdy_fd)
    dCdzstar_fd = 0.5*(dCdx_fd + im*dCdy_fd)

    # auto diff
    dchebt! * C2
    dCdz, dCdzstar = extract_derivative(Nothing,C2)

    @test isapprox(norm(dCdz - dCdz_fd),0.0,atol=TOL)
    @test isapprox(norm(dCdzstar - dCdzstar_fd),0.0,atol=TOL)

    # diff wrt strength
    newblobs = Vortex.dualize_strength(blobs,i,Nothing);
    C2 = zeros(typeof(ComplexRealDual()),N)
    induce_velocity!(C2,p,newblobs,0.0)
    dchebt! * C2
    dCdΓ = extract_derivative(Nothing,C2)

    @test isapprox(norm(dCdΓ - dCdΓ_fd),0.0,atol=TOL)

    # Now with enforce_no_flow_through
    newblobs = Vortex.dualize_position(blobs,i,Nothing)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    dCdz, dCdzstar = extract_derivative(Nothing,pdual.C)

    @test isapprox(norm(dCdz - dCdz_fd),0.0,atol=TOL)
    @test isapprox(norm(dCdzstar - dCdzstar_fd),0.0,atol=TOL)

    n = rand(0:N-1)
    @test isapprox(p.A[n],value(pdual.A[n]),atol=BIGEPS)

    # note that we need to wrap A in complex to ensure it gets dispatched
    # to the correct extract_derivative.
    dAdz, dAdzstar = extract_derivative(Nothing,complex(pdual.A[n]))
    @test dAdz == -0.5im*(dCdz[n+1] - conj(dCdzstar[n+1]))
    @test dAdzstar == conj(dAdz)

    @test extract_derivative(Nothing,pdual.Γ) == 0.0

    # with dualized strength
    newblobs = Vortex.dualize_strength(blobs,i,Nothing);
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    dCdΓ = extract_derivative(Nothing,pdual.C)

    @test isapprox(norm(dCdΓ - dCdΓ_fd),0.0,atol=TOL)

    n = rand(0:N-1)
    @test isapprox(p.A[n],value(pdual.A[n]),atol=BIGEPS)

    @test extract_derivative(Nothing,pdual.Γ) == -1.0
    @test extract_derivative(Nothing,pdual.A[n]) == imag(dCdΓ[n+1])
end

@testset "Induced velocity with bc enforced" begin

    # Now apply the full differentiation to evaluate sensitivity of induced velocity
    # First, by finite difference


    w_fd =  induce_velocity(z,(p,blobs),0.0)
    wx⁺_fd = induce_velocity(z,(px⁺,blobsx⁺),0.0)
    wy⁺_fd = induce_velocity(z,(py⁺,blobsy⁺),0.0)
    wΓ⁺_fd = induce_velocity(z,(pΓ⁺,blobsΓ⁺),0.0)

    dwdx_fd = (wx⁺_fd - w_fd)/dz[i]
    dwdy_fd = (wy⁺_fd - w_fd)/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)
    dwdΓ_fd = (wΓ⁺_fd - w_fd)/dΓ[i]

    # now autodiff
    newblobs = Vortex.dualize_position(blobs,i,Nothing);
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    w = induce_velocity(z,(pdual,newblobs),0.0)

    dwdz, dwdzstar = extract_derivative(Nothing,w)

    @test isapprox(abs(dwdz-dwdz_fd),0.0,atol=TOL)
    @test isapprox(abs(dwdzstar-dwdzstar_fd),0.0,atol=TOL)

    newblobs = Vortex.dualize_strength(blobs,i,Nothing);
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    w = induce_velocity(z,(pdual,newblobs),0.0)
    dwdΓ = extract_derivative(Nothing,w)

    @test isapprox(abs(dwdΓ-dwdΓ_fd),0.0,atol=TOL)
end

@testset "Self-induced velocity with bcs enforced" begin

    # Now check the self-induced velocity
    wself_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wself_fd,blobs, 0.0)
    induce_velocity!(wself_fd, blobs, p, 0.0)

    wselfx⁺_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wselfx⁺_fd,blobsx⁺, 0.0)
    induce_velocity!(wselfx⁺_fd, blobsx⁺, px⁺, 0.0)

    wselfy⁺_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wselfy⁺_fd,blobsy⁺, 0.0)
    induce_velocity!(wselfy⁺_fd, blobsy⁺, py⁺, 0.0)

    wselfΓ⁺_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wselfΓ⁺_fd,blobsΓ⁺, 0.0)
    induce_velocity!(wselfΓ⁺_fd, blobsΓ⁺, pΓ⁺, 0.0)

    dwdx_fd = (wselfx⁺_fd - wself_fd)/dz[i]
    dwdy_fd = (wselfy⁺_fd - wself_fd)/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)
    dwdΓ_fd = (wselfΓ⁺_fd - wself_fd)/dΓ[i]

    newblobs = Vortex.dualize_position(blobs,i,Nothing)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    wself = zeros(typeof(ComplexComplexDual()),length(newblobs))
    self_induce_velocity!(wself,newblobs, 0.0)
    induce_velocity!(wself, newblobs, pdual, 0.0)
    dwdz, dwdzstar  = extract_derivative(Nothing,wself)

    @test isapprox(norm(dwdz-dwdz_fd),0.0,atol=TOL)
    @test isapprox(norm(dwdzstar-dwdzstar_fd),0.0,atol=TOL)

    newblobs = Vortex.dualize_strength(blobs,i,Nothing)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    wself = zeros(typeof(ComplexRealDual()),length(newblobs))
    self_induce_velocity!(wself,newblobs, 0.0)
    induce_velocity!(wself, newblobs, pdual, 0.0)
    dwdΓ = extract_derivative(Nothing,wself)

    @test isapprox(norm(dwdΓ-dwdΓ_fd),0.0,atol=TOL)
end

@testset "Acceleration induced on plate" begin

    targvel = fill(ċ, length(p))

    wself_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wself_fd,blobs, 0.0)
    induce_velocity!(wself_fd, blobs, p, 0.0)
    Ċ_fd = zero(targvel)
    Plates.induce_acc!(Ċ_fd, p.zs, targvel, blobs, wself_fd)

    wselfx⁺_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wselfx⁺_fd,blobsx⁺, 0.0)
    induce_velocity!(wselfx⁺_fd, blobsx⁺, px⁺, 0.0)
    Ċx⁺_fd = zero(targvel)
    Plates.induce_acc!(Ċx⁺_fd, px⁺.zs, targvel, blobsx⁺, wselfx⁺_fd)

    wselfy⁺_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wselfy⁺_fd,blobsy⁺, 0.0)
    induce_velocity!(wselfy⁺_fd, blobsy⁺, py⁺, 0.0)
    Ċy⁺_fd = zero(targvel)
    Plates.induce_acc!(Ċy⁺_fd, py⁺.zs, targvel, blobsy⁺, wselfy⁺_fd)

    wselfΓ⁺_fd = zeros(ComplexF64,length(blobs))
    self_induce_velocity!(wselfΓ⁺_fd,blobsΓ⁺, 0.0)
    induce_velocity!(wselfΓ⁺_fd, blobsΓ⁺, pΓ⁺, 0.0)
    ĊΓ⁺_fd = zero(targvel)
    Plates.induce_acc!(ĊΓ⁺_fd, pΓ⁺.zs, targvel, blobsΓ⁺, wselfΓ⁺_fd)

    dĊdx_fd = (Ċx⁺_fd - Ċ_fd)/dz[i]
    dĊdy_fd = (Ċy⁺_fd - Ċ_fd)/dz[i]
    dĊdz_fd = 0.5*(dĊdx_fd - im*dĊdy_fd)
    dĊdzstar_fd = 0.5*(dĊdx_fd + im*dĊdy_fd)
    dĊdΓ_fd = (ĊΓ⁺_fd - Ċ_fd)/dΓ[i]

    newblobs = Vortex.dualize_position(blobs,i,Nothing);
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0);
    wself = zeros(typeof(ComplexComplexDual()),length(newblobs))
    self_induce_velocity!(wself,newblobs, 0.0)
    induce_velocity!(wself, newblobs, pdual, 0.0)

    Ċ = zeros(typeof(ComplexComplexDual()),length(pdual))
    Plates.induce_acc!(Ċ, pdual.zs, targvel, newblobs, wself)
    dĊdz, dĊdzstar = extract_derivative(Nothing,Ċ)

    @test isapprox(norm(dĊdz-dĊdz_fd),0.0,atol=BIGTOL)
    @test isapprox(norm(dĊdzstar-dĊdzstar_fd),0.0,atol=BIGTOL)

    newblobs = Vortex.dualize_strength(blobs,i,Nothing);
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0);
    wself = zeros(typeof(ComplexRealDual()),length(newblobs))
    self_induce_velocity!(wself,newblobs, 0.0)
    induce_velocity!(wself, newblobs, pdual, 0.0)

    Ċ = zeros(typeof(ComplexRealDual()),length(pdual))
    Plates.induce_acc!(Ċ, pdual.zs, targvel, newblobs, wself)
    dĊdΓ = extract_derivative(Nothing,Ċ)

    @test isapprox(norm(dĊdΓ-dĊdΓ_fd),0.0,atol=BIGTOL)


  end

  Δt = 0.01
  lesp = 0.2
  tesp = 0.0
  z₊ = rand()
  z₋ = rand()

  @testset "Pressure" begin


    press_fd = Plates.surface_pressure_inst(p,motion,blobs,(z₊,z₋),0.0,Δt,lesp,tesp)
    pressx⁺_fd = Plates.surface_pressure_inst(px⁺,motion,blobsx⁺,(z₊,z₋),0.0,Δt,lesp,tesp)
    pressy⁺_fd = Plates.surface_pressure_inst(py⁺,motion,blobsy⁺,(z₊,z₋),0.0,Δt,lesp,tesp)
    pressΓ⁺_fd = Plates.surface_pressure_inst(pΓ⁺,motion,blobsΓ⁺,(z₊,z₋),0.0,Δt,lesp,tesp)

    dpdx_fd = (pressx⁺_fd - press_fd)/dz[i]
    dpdy_fd = (pressy⁺_fd - press_fd)/dz[i]
    dpdz_fd = 0.5*(dpdx_fd - im*dpdy_fd)
    dpdzstar_fd = 0.5*(dpdx_fd + im*dpdy_fd)
    dpdΓ_fd = (pressΓ⁺_fd - press_fd)/dΓ[i]

    newblobs = Vortex.dualize_position(blobs,i,Nothing)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    press = Plates.surface_pressure_inst(pdual,motion,newblobs,(z₊,z₋),0.0,Δt,lesp,tesp)
    dpdz,dpdzstar = extract_derivative(Nothing,complex(press))

    @test isapprox(safenorm(dpdz-dpdz_fd)/safenorm(dpdz),0.0,atol=BIGTOL)
    @test isapprox(safenorm(dpdzstar-dpdzstar_fd)/safenorm(dpdzstar),0.0,atol=BIGTOL)

    newblobs = Vortex.dualize_strength(blobs,i,Nothing)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    press = Plates.surface_pressure_inst(pdual,motion,newblobs,(z₊,z₋),0.0,Δt,lesp,tesp)
    dpdΓ = extract_derivative(Nothing,press)

    @test isapprox(safenorm(dpdΓ-dpdΓ_fd)/safenorm(dpdΓ),0.0,atol=BIGTOL)


  end

end
