using LinearAlgebra

import PotentialFlow.Utils: derivative, extract_derivative, value,
          Dual,ComplexComplexDual,ComplexRealDual

const BIGEPS = 1000*eps(1.0)
const TOL=5e-6

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

  @testset "Basic operations with duals" begin

    pos = rand(ComplexF64,5)
    str = rand(length(pos))

    σ = 1e-2
    blobs = Vortex.Blob.(pos,str,σ)
    z = rand(ComplexF64)

    i = 3


    dualpos = ComplexComplexDual{Nothing}(Elements.position(blobs)[i],one(z),zero(z))

    newblob = Vortex.Blob(dualpos,Elements.circulation(blobs[i]),σ)
    @test value(Elements.position(newblob)) == Elements.position(blobs)[i]

    dwdz,dwdzstar  = extract_derivative(Nothing,induce_velocity(z,newblob,0.0))

    w2 = Elements.circulation(newblob)*PotentialFlow.Blobs.blob_kernel(z - Elements.position(newblob),σ)
    dwdz2, dwdzstar2 = extract_derivative(Nothing,w2)

    @test isapprox(abs(dwdz-dwdz2),0.0,atol=BIGEPS)
    @test isapprox(abs(dwdzstar-dwdzstar2),0.0,atol=BIGEPS)

    # Finite difference approximations
    eps = 1e-10
    dz = zeros(ComplexF64,length(blobs))
    dz[i] = eps
    blobsx⁺ = Vortex.Blob.(Elements.position(blobs).+dz,Elements.circulation.(blobs),σ)
    dwdx_fd = (induce_velocity(z,blobsx⁺,0.0) - induce_velocity(z,blobs,0.0))/dz[i]

    blobsy⁺ = Vortex.Blob.(Elements.position(blobs).+im*dz,Elements.circulation.(blobs),σ)
    dwdy_fd =(induce_velocity(z,blobsy⁺,0.0) - induce_velocity(z,blobs,0.0))/dz[i]

    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)

    dΓ = zeros(Float64,length(blobs))
    dΓ[i] = eps
    blobsΓ⁺ = Vortex.Blob.(Elements.position(blobs),Elements.circulation.(blobs).+dΓ,σ)
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


    N = 7
    C  = zeros(ComplexF64, N)
    dchebt! = Plates.Chebyshev.plan_transform!(C)

    p = PotentialFlow.Plate(N,2.0,complex(0),0.0)
    motion = PotentialFlow.RigidBodyMotion(complex(0),0.0)

    C  = zeros(ComplexF64, N)
    induce_velocity!(C,p,blobs,0.0)

    #blobsx⁺ = Vortex.Blob.(Elements.position(blobs).+dz,Elements.circulation.(blobs),σ);
    Cx⁺  = zeros(ComplexF64, N)
    induce_velocity!(Cx⁺,p,blobsx⁺,0.0)
    dwdx_fd = (Cx⁺ - C)/dz[i]

    #blobsy⁺ = Vortex.Blob.(Elements.position(blobs).+im*dz,Elements.circulation.(blobs),σ);
    Cy⁺  = zeros(ComplexF64, N)
    induce_velocity!(Cy⁺,p,blobsy⁺,0.0)
    dwdy_fd = (Cy⁺ - C)/dz[i]

    #blobsΓ⁺ = Vortex.Blob.(Elements.position(blobs),Elements.circulation.(blobs).+dΓ,σ);
    CΓ⁺  = zeros(ComplexF64, N)
    induce_velocity!(CΓ⁺,p,blobsΓ⁺,0.0)

    dwdz_fd = 0.5*(dwdx_fd - im*dwdy_fd)
    dwdzstar_fd = 0.5*(dwdx_fd + im*dwdy_fd)

    newblobs = Vortex.dualize_position(blobs,i,Nothing)
    C2 = zeros(typeof(complex(Dual{Nothing}(0.0,0.0,0.0))),N)
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
    C2 = zeros(typeof(complex(Dual{Nothing}(0.0,0.0))),N)
    induce_velocity!(C2,p,newblobs,0.0)
    dchebt! * C2
    dCdΓ = extract_derivative(Nothing,C2)

    @test isapprox(norm(dCdΓ - dCdΓ_fd),0.0,atol=TOL)


  end

end
