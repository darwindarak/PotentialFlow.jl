using ForwardDiff

@testset "Complex Automatic Differentiation" begin

  z = rand(ComplexF64)

  dz1, dzstar1 = ForwardDiff.derivative(z -> log(sqrt(z)),z)

  dzex = 1/(2*z)

  @test isapprox(abs(dz1-dzex),0,atol=10*eps(1.0))
  @test isapprox(abs(dzstar1),0,atol=10*eps(1.0))

  dz1, dzstar1 = ForwardDiff.derivative(z -> conj(z),z)
  @test dzstar1 == one(z) && dz1 == zero(z)

  f(z) = 0.5im/(π*conj(z))

  dz1, dzstar1 = ForwardDiff.derivative(f,z)
  dzstarex = -0.5im/(pi*conj(z)^2)

  @test isapprox(abs(dz1),0,atol=10*eps(1.0))
  @test isapprox(abs(dzstar1-dzstarex),0,atol=10*eps(1.0))

  z0 = rand(ComplexF64)
  dz1, dzstar1 = ForwardDiff.derivative(z -> f(z-z0),z)

  dzstarex = -0.5im/(pi*conj(z-z0)^2)

  @test isapprox(abs(dz1),0,atol=10*eps(1.0))
  @test isapprox(abs(dzstar1-dzstarex),0,atol=10*eps(1.0))

  dz1, dzstar1 = ForwardDiff.derivative(abs,z)

  dzex = 0.5*conj(z)/abs(z)
  dzstarex = 0.5*z/abs(z)
  @test isapprox(abs(dz1-dzex),0,atol=10*eps(1.0))
  @test isapprox(abs(dzstar1-dzstarex),0,atol=10*eps(1.0))

  dz1, dzstar1 = ForwardDiff.derivative(z -> z + sqrt(log(z)) - 1/z^2,z)

  dzex = 1 + 1/(2*sqrt(log(z))*z) + 2/z^3
  dzstarex = 0.0

  @test isapprox(abs(dz1-dzex),0,atol=10*eps(1.0))
  @test isapprox(abs(dzstar1-dzstarex),0,atol=10*eps(1.0))

  dz1, dzstar1 = ForwardDiff.derivative(z -> real(z - √(z - 1)*√(z + 1)),z)

  dzex = 0.5*(1 - 0.5*sqrt((z+1)/(z-1)) - 0.5*sqrt((z-1)/(z+1)))
  dzstarex = 0.5*conj(1 - 0.5*sqrt((z+1)/(z-1)) - 0.5*sqrt((z-1)/(z+1)))

  @test isapprox(abs(dz1-dzex),0,atol=10*eps(1.0))
  @test isapprox(abs(dzstar1-dzstarex),0,atol=10*eps(1.0))

end
