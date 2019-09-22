using Test, LinearAlgebra, Mire


@testset "QG Inertial modes" begin
    a,b,c=1.,1.,1.
    Ω = [0,0,1]
    N = 3
    B,A, uj = assemblehd(N, a, b, c, Ω)
    esol = eigen(inv(Matrix(B))*A)
    ω = esol.values

    # Zhang et al., J. Fluid Mech. (2001), vol. 437, pp. 103–119. eq (4.5)
    zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im
    @test any(esol.values.≈zhang(1,1))
    @test any(esol.values.≈zhang(2,1))
end
