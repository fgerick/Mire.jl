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

@testset "QG Malkus modes" begin
    a,b,c=1.,1.,1.
    Le=1e-6
    Ω = [0,0,1/Le]
    N = 3
    b0 = [-y,x,0]
    B,A, uj = assemblemhd(N, a, b, c, Ω,b0)
    esol = eigen(inv(Matrix(B))*A)
    ω = esol.values

    # Zhang et al., J. Fluid Mech. (2001), vol. 437, pp. 103–119, eq. (4.5)
    zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im

    # Malkus J. Fluid Mech. (1967), vol. 28, pp. 793-802, eq. (2.28)
    slow(m, N, Le, λ = imag(zhang(m,N))) = im*λ/2Le*(1 - √(1+4Le^2*m*(m-λ)/λ^2))
    fast(m, N, Le, λ = imag(zhang(m,N))) = im*λ/2Le*(1 + √(1+4Le^2*m*(m-λ)/λ^2))

    @test any(isapprox.(esol.values,slow(1,1,Le),atol=1e-10))
    @test any(isapprox.(esol.values,slow(2,1,Le),atol=1e-10))
    @test any(isapprox.(esol.values,fast(1,1,Le),atol=1e-9))
    @test any(isapprox.(esol.values,fast(2,1,Le),atol=1e-9))
end

@testset "QG from scalar vs QG 3D projection" begin
    N,a,b,c,Le = 5,1.25,0.8,1.0,1e-5
    Ω = 1/Le
    A0 = (x^0*y^0 + x)/3
    b0 = (Mire.qg_vel(0,0,a,b,c) + Mire.qg_vel(1,0,a,b,c))/3

    LHS,RHS,vsqg=Mire.assemblemhd_quag(N,a,b,c,Ω,A0)
    LHS2,RHS2,vsqg2 = Mire.assemblemhd_qg(N,a,b,c,[0,0,Ω],b0)

    s1 = eigen(Matrix(RHS),Matrix(LHS))
    s2 = eigen(Matrix(RHS2),Matrix(LHS2))

    @test sort(abs.(s1.values),rev=true) ≈ sort(abs.(s2.values),rev=true)

end
