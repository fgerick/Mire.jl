

@testset "Columnar 3D inertial modes" begin
    # a,b,c=1.,1.,1.
    Ω = [0,0,1.0]
    N = 3
    prob = HDProblem(N,Sphere(), Ω, LebovitzBasis)
    assemble!(prob)
    A,B = prob.RHS,prob.LHS
    esol = eigen(inv(Matrix(B))*A)
    ω = esol.values

    # Zhang et al., J. Fluid Mech. (2001), vol. 437, pp. 103–119. eq (4.5)
    zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im
    @test any(esol.values.≈zhang(1,1))
    @test any(esol.values.≈zhang(2,1))
end

@testset "Columnar 3D Malkus modes" begin
    Le = 1e-6
    Ω = [0,0,1/Le]
    N = 3
    b0 = [-y,x,0]
    prob = MHDProblem(N, Sphere(), Ω, b0, LebovitzBasis, ConductingMFBasis)
    assemble!(prob)
    A,B = prob.RHS,prob.LHS
    esol = eigen(inv(Matrix(B))*A)
    ω = esol.values

    # Zhang et al., J. Fluid Mech. (2001), vol. 437, pp. 103–119, eq. (4.5)
    zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im
    # note: for larger m,N this approximate formula is inaccurate, and the comparison will fail
    # For the exact values of the inertial mode frequencies the roots of the univariate polynomials
    # given in (2.15/2.16) should be used as a comparison to this code.


    # Malkus J. Fluid Mech. (1967), vol. 28, pp. 793-802, eq. (2.28)
    slow(m, N, Le, λ = imag(zhang(m,N))) = im*λ/2Le*(1 - √(1+4Le^2*m*(m-λ)/λ^2))
    fast(m, N, Le, λ = imag(zhang(m,N))) = im*λ/2Le*(1 + √(1+4Le^2*m*(m-λ)/λ^2))

    @test any(isapprox.(esol.values,slow(1,1,Le),atol=1e-10))
    @test any(isapprox.(esol.values,slow(2,1,Le),atol=1e-10))
    @test any(isapprox.(esol.values,fast(1,1,Le),atol=1e-9))
    @test any(isapprox.(esol.values,fast(2,1,Le),atol=1e-9))
end

@testset "Spheroidal QG inertial modes" begin
    a,b = 1.0, 1.1
    Ω = [0,0,1.0]
    N = 9
    prob = HDProblem(N, Ellipsoid(a,a,b), Ω, QGBasis)
    assemble!(prob)
    A,B = prob.RHS,prob.LHS
    esol = eigen(inv(Matrix(B))*A)
    ω = esol.values

    # Maffei et al. 2017, Proc. R. Soc. A 473: 20170181. http://dx.doi.org/10.1098/rspa.2017.0181 eq. (4.11)

    maffei(m,N,b)= -m/(N*(2N + 2m + 1) + m/2 + m^2*b^2/6)*im
    for m in 1:3, N in 1:3
        @test any(esol.values.≈maffei(m,N,b)) 
    end
end

@testset "λ₁₂ ellipsoidal inertial mode" begin
    a,b,c = 1.2,1.1,0.7
    Ω = [0,0,1.0]
    N = 5
    prob = HDProblem(N, Ellipsoid(a,b,c), Ω, LebovitzBasis)
    assemble!(prob)
    A,B = prob.RHS,prob.LHS
    esol = eigen(inv(Matrix(B))*A)
    ω = esol.values

    # Vantieghem 2014, Proc. R. Soc. A 470: 20140093. http://dx.doi.org/10.1098/rspa.2014.0093 eq. (3.21)
    vantieghem(a,b,c)= 2a*b/sqrt((a^2+c^2)*(b^2+c^2))*im

    @test any(esol.values.≈vantieghem(a,b,c))
    @test any(esol.values.≈-vantieghem(a,b,c))  
end

@testset "QG from scalar vs QG 3D projection" begin
 
    N,a,b,c,Le = 5,1.25,0.8,0.9,1e-5  
    Ω = [0,0,1/Le]
    V = Ellipsoid(a,b,c)
    b0 = (Mire.uqg(0,0,V) + Mire.uqg(1,0,V))/3

    prob = MHDProblem(N, V, Ω, b0, QGBasis, QGBasis)
    assemble!(prob)
    RHS,LHS = prob.RHS,prob.LHS 

    Ω = 1/Le
    A0 = (x^0*y^0 + x)/3

    LHSquag,RHSquag,vsqg=Mire.Quagmire.assemblemhd_quag(N,a,b,c,Ω,A0)

    s1 = eigen(Matrix(RHS),Matrix(LHS))
    s2 = eigen(Matrix(RHSquag),Matrix(LHSquag))

    @test sort(abs.(s1.values),rev=true) ≈ sort(abs.(s2.values),rev=true)

end
