using LinearAlgebra

import PotentialFlow.Utils: derivative, extract_derivative, value, partials,
          Dual,ComplexDual, ComplexGradientConfig,
          dz_partials, gradient, chunk_mode_gradient, chunk_mode_jacobian,
          complexderivative

import PotentialFlow.Elements: gradient_position, gradient_strength,
          jacobian_position, jacobian_strength, jacobian_param

const DELTA=1e-6
const BIGEPS = 1e-8
const TOL=5e-6
const BIGTOL=1e-4
const BIGGESTTOL=1e-3

safenorm(a) = norm(filter(x -> ~isnan(x),a))

@testset "Complex Automatic Differentiation" begin

  @testset "Basic derivatives" begin

    z = rand(ComplexF64)

    dz1, dzstar1 = complexderivative(z -> log(sqrt(z)),z)

    dzex = 1/(2*z)

    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1),0,atol=BIGEPS)

    dz1, dzstar1 = complexderivative(z -> conj(z),z)
    @test dzstar1 == one(z) && dz1 == zero(z)

    f(z) = 0.5im/(π*conj(z))

    dz1, dzstar1 = complexderivative(f,z)
    dzstarex = -0.5im/(pi*conj(z)^2)

    @test isapprox(abs(dz1),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    z0 = rand(ComplexF64)
    dz1, dzstar1 = complexderivative(z -> f(z-z0),z)

    dzstarex = -0.5im/(pi*conj(z-z0)^2)

    @test isapprox(abs(dz1),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    dz1, dzstar1 = complexderivative(abs,z)

    dzex = 0.5*conj(z)/abs(z)
    dzstarex = 0.5*z/abs(z)
    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    dz1, dzstar1 = complexderivative(z -> z + sqrt(log(z)) - 1/z^2,z)

    dzex = 1 + 1/(2*sqrt(log(z))*z) + 2/z^3
    dzstarex = 0.0

    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

    dz1, dzstar1 = complexderivative(z -> real(z - √(z - 1)*√(z + 1)),z)

    dzex = 0.5*(1 - 0.5*sqrt((z+1)/(z-1)) - 0.5*sqrt((z-1)/(z+1)))
    dzstarex = 0.5*conj(1 - 0.5*sqrt((z+1)/(z-1)) - 0.5*sqrt((z-1)/(z+1)))

    @test isapprox(abs(dz1-dzex),0,atol=BIGEPS)
    @test isapprox(abs(dzstar1-dzstarex),0,atol=BIGEPS)

  end

  @testset "Nested derivatives" begin

    z0 = rand(ComplexF64)
    y0 = rand(ComplexF64)
    dyz, dystarz, dyzstar, dystarzstar = complexderivative(y -> complexderivative(z -> 1/(z-y),z0),y0)

    dyzex = -2/(z0-y0)^3
    @test isapprox(abs(dyz-dyzex),0,atol=BIGEPS)
    @test isapprox(abs(dystarz),0,atol=BIGEPS)
    @test isapprox(abs(dyzstar),0,atol=BIGEPS)
    @test isapprox(abs(dystarzstar),0,atol=BIGEPS)

    z0 = rand(ComplexF64)
    y0 = rand(ComplexF64)
    dyz, dystarz, dyzstar, dystarzstar = complexderivative(y -> complexderivative(z -> 1/(z-conj(y)),z0),y0)

    dystarzex = -2/(z0-conj(y0))^3
    @test isapprox(abs(dyz),0,atol=BIGEPS)
    @test isapprox(abs(dystarz-dystarzex),0,atol=BIGEPS)
    @test isapprox(abs(dyzstar),0,atol=BIGEPS)
    @test isapprox(abs(dystarzstar),0,atol=BIGEPS)

    z0 = rand(ComplexF64)
    y0 = rand(ComplexF64)
    dyz, dystarz, dyzstar, dystarzstar = complexderivative(y -> complexderivative(z -> 1/(conj(z)-y),z0),y0)

    dyzstarex = -2/(conj(z0)-y0)^3
    @test isapprox(abs(dyz),0,atol=BIGEPS)
    @test isapprox(abs(dystarz),0,atol=BIGEPS)
    @test isapprox(abs(dyzstar-dyzstarex),0,atol=BIGEPS)
    @test isapprox(abs(dystarzstar),0,atol=BIGEPS)

    z0 = rand(ComplexF64)
    y0 = rand(ComplexF64)
    dyz, dystarz, dyzstar, dystarzstar = complexderivative(y -> complexderivative(z -> 1/conj(z-y),z0),y0)

    dystarzstarex = -2/conj(z0-y0)^3
    @test isapprox(abs(dyz),0,atol=BIGEPS)
    @test isapprox(abs(dystarz),0,atol=BIGEPS)
    @test isapprox(abs(dyzstar),0,atol=BIGEPS)
    @test isapprox(abs(dystarzstar-dystarzstarex),0,atol=BIGEPS)

  end

  @testset "Chunk mode derivatives" begin

    n = 20
    z = rand(ComplexF64,n)
    fcn(z) = sum(log.(sqrt.(z)))

    cfg = ComplexGradientConfig(fcn,z)
    dz, dzstar = gradient(fcn,z,cfg)

    ichk = rand(1:n)
    @test dz[ichk] ≈ 1.0/(2*z[ichk])

    N = rand(1:n)
    cfg2 = ComplexGradientConfig(fcn,z,PotentialFlow.Utils.Chunk{N}())
    dz_chunk, dzstar_chunk = gradient(fcn,z,cfg2)
    @test dz == dz_chunk
    @test dzstar == dzstar_chunk

    n = 20
    z = rand(ComplexF64,n)
    fcn2(z) = log.(sqrt.(z))

    cfg = ComplexGradientConfig(fcn2,z)
    dz, dzstar = PotentialFlow.Utils.jacobian(fcn2,z,cfg)

    ichk = rand(1:n)
    @test dz[ichk,ichk] ≈ 1.0/(2*z[ichk])

    N = rand(1:n)
    cfg2 = ComplexGradientConfig(fcn2,z,PotentialFlow.Utils.Chunk{N}())

    dz_chunk, dzstar_chunk = PotentialFlow.Utils.jacobian(fcn2,z,cfg2)

    @test dz_chunk == dz
    @test dzstar_chunk == dzstar

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


    dualpos = one(ComplexDual{Nothing},Elements.position(blobs)[i])

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

    compute_velocity(v) = induce_velocity(z,v,0.0)

    # Finite difference approximations
    w_fd = compute_velocity(blobs)
    wx⁺_fd = compute_velocity(blobsx⁺)
    wy⁺_fd = compute_velocity(blobsy⁺)
    wΓ⁺_fd = compute_velocity(blobsΓ⁺)
    dwdx_fd = (wx⁺_fd - w_fd)/dz[i]
    dwdy_fd = (wy⁺_fd - w_fd)/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)
    dwdΓ_fd = (wΓ⁺_fd - w_fd)/dΓ[i]

    #=
    dwdx_fd = (induce_velocity(z,blobsx⁺,0.0) - induce_velocity(z,blobs,0.0))/dz[i]
    dwdy_fd =(induce_velocity(z,blobsy⁺,0.0) - induce_velocity(z,blobs,0.0))/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)
    dwdΓ_fd = (induce_velocity(z,blobsΓ⁺,0.0) - induce_velocity(z,blobs,0.0))/dΓ[i]
    =#

    # Auto differentation
    cfg = ComplexGradientConfig(z -> (),blobs)
    newblobs = Vortex.seed_position(blobs,cfg)
    dwdz2, dwdzstar2 = dz_partials(induce_velocity(z,newblobs,0.0))

    # chunk mode
    nchunk = 3
    cfg2 = PotentialFlow.Utils.ComplexGradientConfig(z -> (),blobs,PotentialFlow.Utils.Chunk{nchunk}())

    # Auto diff with the API
    dwdz, dwdzstar = gradient_position(compute_velocity,blobs)
    @test dwdz == dwdz2 && dwdzstar == dwdzstar2

    @test_skip abs(dwdz[i]-dwdz_fd) ≈ 0.0 atol=TOL
    @test_skip abs(dwdzstar[i]-dwdzstar_fd) ≈ 0.0 atol=TOL

    dwdz_chunk, dwdzstar_chunk = gradient_position(compute_velocity,blobs,cfg2)
    @test dwdz == dwdz_chunk && dwdzstar == dwdzstar_chunk

    newblobs = Vortex.seed_strength(blobs,cfg)

    @test sum(value.(Vortex.circulation.(newblobs))) -
              value(Vortex.circulation(newblobs)) == 0

    dwdz_tmp, dwdzstar_tmp = dz_partials(induce_velocity(z,newblobs,0.0))
    dwdΓ2 = dwdz_tmp+dwdzstar_tmp

    # Auto diff with the API
    dwdΓ = gradient_strength(compute_velocity,blobs)
    @test dwdΓ == dwdΓ2

    @test_skip abs(dwdΓ[i]-dwdΓ_fd) ≈ 0.0 atol=TOL

    dwdΓ_chunk = gradient_strength(compute_velocity,blobs,cfg2)

    @test dwdΓ  == dwdΓ_chunk

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

    function plate_velocity(v)
      w = zeros(Complex{Elements.property_type(eltype(v))},N)
      induce_velocity!(w,p,v,0.0)
      return w
    end

    C = plate_velocity(blobs)
    Cx⁺ = plate_velocity(blobsx⁺)
    Cy⁺ = plate_velocity(blobsy⁺)
    CΓ⁺ = plate_velocity(blobsΓ⁺)


    dwdx_fd = (Cx⁺ - C)/dz[i]
    dwdy_fd = (Cy⁺ - C)/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)

    # Auto diff
    dwdz, dwdzstar = jacobian_position(plate_velocity,blobs)

    # test that the induced velocities and their derivatives match
    #@test isapprox(norm(value.(C2) - C),0.0,atol=BIGEPS)
    @test_skip norm(dwdz[:,i] - dwdz_fd) ≈ 0.0 atol=TOL
    @test_skip norm(dwdzstar[:,i] - dwdzstar_fd) ≈ 0.0 atol=TOL

    dchebt! = Plates.Chebyshev.plan_transform!(Plates._dct_data(Float64,N))


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
    cfg = ComplexGradientConfig(z -> (),blobs)
    newblobs = Vortex.seed_position(blobs,cfg)
    C2 = zeros(eltype(cfg.duals),N)
    induce_velocity!(C2,p,newblobs,0.0)
    dchebt! = Plates.Chebyshev.plan_transform!(Plates._dct_data(Dual{Nothing,Float64,2*nblob},N))
    dchebt! * C2
    dCdz, dCdzstar = dz_partials(C2,i)

    @test_skip norm(dCdz - dCdz_fd) ≈ 0.0 atol=TOL
    @test_skip norm(dCdzstar - dCdzstar_fd) ≈ 0.0 atol=TOL

    # diff wrt strength
    newblobs = Vortex.seed_strength(blobs,cfg)
    C2 = zeros(eltype(cfg.duals),N)
    induce_velocity!(C2,p,newblobs,0.0)
    dchebt! * C2
    dCdz_tmp, dCdzstar_tmp = dz_partials(C2,i)
    dCdΓ = dCdz_tmp+dCdzstar_tmp

    @test_skip norm(dCdΓ - dCdΓ_fd) ≈ 0.0 atol=TOL

    # Now with enforce_no_flow_through
    newblobs = Vortex.seed_position(blobs,cfg)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    dCdz, dCdzstar = dz_partials(pdual.C,i)

    @test_skip norm(dCdz - dCdz_fd) ≈ 0.0 atol=TOL
    @test_skip norm(dCdzstar - dCdzstar_fd) ≈ 0.0 atol=TOL

    n = rand(0:N-1)
    @test p.A[n] ≈ value(pdual.A[n]) atol=BIGEPS

    # note that we need to wrap A in complex to ensure it gets dispatched
    # to the correct extract_derivative.
    dAdz, dAdzstar = dz_partials(complex(pdual.A[n]),i)
    @test dAdz == -0.5im*(dCdz[n+1] - conj(dCdzstar[n+1]))
    @test dAdzstar == conj(dAdz)

    dΓdz, dΓdzstar = dz_partials(complex(pdual.Γ))
    @test all(dΓdz .== complex(0.0))

    # with dualized strength
    newblobs = Vortex.seed_strength(blobs,cfg)
    pdual = PotentialFlow.Plate{Elements.property_type(eltype(newblobs))}(N,L,c,α)
    Plates.enforce_no_flow_through!(pdual, motion, newblobs, 0.0)
    dCdz_tmp, dCdzstar_tmp = dz_partials(pdual.C,i)
    dCdΓ = dCdz_tmp+dCdzstar_tmp

    @test_skip norm(dCdΓ - dCdΓ_fd) ≈ 0.0 atol=TOL

    n = rand(0:N-1)
    @test p.A[n] ≈ value(pdual.A[n]) atol=BIGEPS

    dΓdz, dΓdzstar = dz_partials(complex(pdual.Γ))
    dΓdΓ = dΓdz+dΓdzstar
    @test all(dΓdΓ .== complex(-1.0))

    @test pdual.A[n] == imag(pdual.C[n+1])
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
    function f(v)
      ptmp = PotentialFlow.Plate{Elements.property_type(eltype(v))}(N,L,c,α)
      Plates.enforce_no_flow_through!(ptmp, motion, v, 0.0)
      out = induce_velocity(z,(ptmp,v),0.0)
      return out
    end

    dwdz, dwdzstar = gradient_position(f,blobs)

    @test_skip abs(dwdz[i]-dwdz_fd) ≈ 0.0 atol=TOL
    @test_skip abs(dwdzstar[i]-dwdzstar_fd) ≈ 0.0 atol=TOL

    dwdΓ = gradient_strength(f,blobs)

    @test_skip abs(dwdΓ[i]-dwdΓ_fd) ≈ 0.0 atol=TOL
end


@testset "Self-induced velocity with bcs enforced" begin

    function self_velocity(v)
      ptmp = PotentialFlow.Plate{Elements.property_type(eltype(v))}(N,L,c,α)
      Plates.enforce_no_flow_through!(ptmp, motion, v, 0.0)
      wself = zeros(Complex{Elements.property_type(eltype(v))},length(v))
      self_induce_velocity!(wself,v, 0.0)
      induce_velocity!(wself, v, ptmp, 0.0)
      return wself
    end

    # chunk mode
    nchunk = 3
    cfg2 = PotentialFlow.Utils.ComplexGradientConfig(z -> (),blobs,PotentialFlow.Utils.Chunk{nchunk}())


    wself_fd = self_velocity(blobs)
    wselfx⁺_fd = self_velocity(blobsx⁺)
    wselfy⁺_fd = self_velocity(blobsy⁺)
    wselfΓ⁺_fd = self_velocity(blobsΓ⁺)
    dwdx_fd = (wselfx⁺_fd - wself_fd)/dz[i]
    dwdy_fd = (wselfy⁺_fd - wself_fd)/dz[i]
    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)
    dwdΓ_fd = (wselfΓ⁺_fd - wself_fd)/dΓ[i]

    dwdz, dwdzstar = jacobian_position(self_velocity,blobs)


    @test_skip norm(dwdz[:,i]-dwdz_fd) ≈ 0.0 atol=TOL
    @test_skip norm(dwdzstar[:,i]-dwdzstar_fd) ≈ 0.0 atol=TOL

    dwdz_chunk, dwdzstar_chunk = jacobian_position(self_velocity,blobs,cfg2)
    @test dwdz == dwdz_chunk
    @test dwdzstar == dwdzstar_chunk


    dwdΓ = jacobian_strength(self_velocity,blobs)

    @test_skip norm(dwdΓ[:,i]-dwdΓ_fd) ≈ 0.0 atol=TOL

    dwdΓ_chunk = jacobian_strength(self_velocity,blobs,cfg2)
    @test dwdΓ == dwdΓ_chunk

end


@testset "Acceleration induced on plate" begin

    targvel = fill(ċ, length(p))

    function compute_Ċ(v)
      ptmp = PotentialFlow.Plate{Elements.property_type(eltype(v))}(N,L,c,α)
      Plates.enforce_no_flow_through!(ptmp, motion, v, 0.0)
      wself = zeros(Complex{Elements.property_type(eltype(v))},length(v))
      self_induce_velocity!(wself,v, 0.0)
      induce_velocity!(wself, v, ptmp, 0.0)
      Ċ = zeros(Complex{Elements.property_type(eltype(v))},length(ptmp))
      Plates.induce_acc!(Ċ, ptmp.zs, targvel, v, wself)
    end

    Ċ_fd = compute_Ċ(blobs)
    Ċx⁺_fd = compute_Ċ(blobsx⁺)
    Ċy⁺_fd = compute_Ċ(blobsy⁺)
    ĊΓ⁺_fd = compute_Ċ(blobsΓ⁺)


    dĊdx_fd = (Ċx⁺_fd - Ċ_fd)/dz[i]
    dĊdy_fd = (Ċy⁺_fd - Ċ_fd)/dz[i]
    dĊdz_fd = 0.5*(dĊdx_fd - im*dĊdy_fd)
    dĊdzstar_fd = 0.5*(dĊdx_fd + im*dĊdy_fd)
    dĊdΓ_fd = (ĊΓ⁺_fd - Ċ_fd)/dΓ[i]


    dĊdz, dĊdzstar = jacobian_position(compute_Ċ,blobs)

    @test_skip norm(dĊdz[:,i]-dĊdz_fd) ≈ 0.0 atol=BIGTOL
    @test_skip norm(dĊdzstar[:,i]-dĊdzstar_fd) ≈ 0.0 atol=BIGTOL


    dĊdΓ = jacobian_strength(compute_Ċ,blobs)

    @test_skip norm(dĊdΓ[:,i]-dĊdΓ_fd) ≈ 0.0 atol=BIGTOL


  end

  Δt = 0.01
  lesp = 0.2
  tesp = 0.0
  z₊ = rand()
  z₋ = rand()

  plesp⁺ = deepcopy(p)
  dlesp = DELTA

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

    function compute_pressure(v)
      ptmp = PotentialFlow.Plate{Elements.property_type(eltype(v))}(N,L,c,α)
      return complex(Plates.surface_pressure_inst(ptmp,motion,v,(z₊,z₋),0.0,Δt,lesp,tesp))
    end

    dpdz, dpdzstar = jacobian_position(compute_pressure,blobs)

    @test safenorm(dpdz[:,i]-dpdz_fd)/safenorm(dpdz_fd) ≈ 0.0 atol=BIGTOL
    @test safenorm(dpdzstar[:,i]-dpdzstar_fd)/safenorm(dpdzstar_fd) ≈ 0.0 atol=BIGTOL


    dpdΓ = real(jacobian_strength(compute_pressure,blobs))

    @test safenorm(dpdΓ[:,i]-dpdΓ_fd)/safenorm(dpdΓ_fd) ≈ 0.0 atol=BIGTOL

    presslesp⁺_fd = Plates.surface_pressure_inst(plesp⁺,motion,blobs,(z₊,z₋),0.0,Δt,lesp+dlesp,tesp)
    dpdlesp_fd = (presslesp⁺_fd - press_fd)/dlesp

    function lesp_to_pressure(v,lesp)
      ptmp = PotentialFlow.Plate{Elements.property_type(eltype(v))}(N,L,c,α)
      press = Plates.surface_pressure_inst(ptmp,motion,v,(z₊,z₋),0.0,Δt,lesp,tesp)
      return complex(press)
    end

    dpdlesp = jacobian_param(lesp_to_pressure,(blobs,lesp))

    @test safenorm(dpdlesp-dpdlesp_fd)/safenorm(dpdlesp) ≈ 0.0 atol=BIGTOL


  end

  ## Larger number of blobs

  nblob = 20
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

  @testset "Acceleration induced on plate" begin

      targvel = fill(ċ, length(p))

      function compute_Ċ(v)
        ptmp = PotentialFlow.Plate{Elements.property_type(eltype(v))}(N,L,c,α)
        Plates.enforce_no_flow_through!(ptmp, motion, v, 0.0)
        wself = zeros(Complex{Elements.property_type(eltype(v))},length(v))
        self_induce_velocity!(wself,v, 0.0)
        induce_velocity!(wself, v, ptmp, 0.0)
        Ċ = zeros(Complex{Elements.property_type(eltype(v))},length(ptmp))
        Plates.induce_acc!(Ċ, ptmp.zs, targvel, v, wself)
      end

      Ċ_fd = compute_Ċ(blobs)
      Ċx⁺_fd = compute_Ċ(blobsx⁺)
      Ċy⁺_fd = compute_Ċ(blobsy⁺)
      ĊΓ⁺_fd = compute_Ċ(blobsΓ⁺)


      dĊdx_fd = (Ċx⁺_fd - Ċ_fd)/dz[i]
      dĊdy_fd = (Ċy⁺_fd - Ċ_fd)/dz[i]
      dĊdz_fd = 0.5*(dĊdx_fd - im*dĊdy_fd)
      dĊdzstar_fd = 0.5*(dĊdx_fd + im*dĊdy_fd)
      dĊdΓ_fd = (ĊΓ⁺_fd - Ċ_fd)/dΓ[i]


      dĊdz, dĊdzstar = jacobian_position(compute_Ċ,blobs)

      @test_skip norm(dĊdz[:,i]-dĊdz_fd) ≈ 0.0 atol=BIGTOL
      @test_skip norm(dĊdzstar[:,i]-dĊdzstar_fd) ≈ 0.0 atol=BIGTOL


      dĊdΓ = jacobian_strength(compute_Ċ,blobs)

      @test_skip norm(dĊdΓ[:,i]-dĊdΓ_fd) ≈ 0.0 atol=BIGTOL


    end

end
