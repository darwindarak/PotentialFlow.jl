@testset "Pressure" begin

function p(z::ComplexF64;Γv=0.0,zv=complex(0.0))
  -Γv^2/(8π^2)*abs2(1/(z-zv)+1/(z+zv)) - Γv^2/(8π^2)*real(1/conj(zv)*(1/(z+zv)-1/(z-zv)))
end

zv = rand(ComplexF64)
Γv = rand(Float64)
ze = rand(ComplexF64)
blobs = Vortex.Blob.([zv,-zv],[Γv,Γv],0.00001)
pts = Vortex.Blob.([zv,-zv],[Γv,Γv])

@test isapprox(p(ze,Γv=Γv,zv=zv)[1],pressure([ze],blobs,0.0)[1],atol=1e-6)
@test isapprox(p(ze,Γv=Γv,zv=zv)[1],pressure([ze],pts,0.0))


end
