@testset "Miscellaneous API" begin


    vs = [[x,y,z],Mire.ex,Mire.ey,Mire.ez]
    bs = [[x,y,z],Mire.ex,Mire.ey,Mire.ez] 
    αs = [1//2,1,0,3, 1,2,0,1]

    @test velocities(vs,αs)[1] == [1//2*x+1,1//2*y,1//2*z+3]
    @test magneticfields(vs,αs)[1] == [x+2//1, y, z + 1//1]

end


@testset "Monomial integrations" begin

    a,b,c = 0.2,0.3,0.4
    cmat = cacheint(3,Ellipsoid(a,b,c))

    @test cmat[1,1,1] ≈ int_monomial_ellipsoid(x^0,a,b,c) ≈ 4/3*a*b*c
    @test cmat[3,5,3] ≈ int_monomial_ellipsoid(x^2*y^4*z^2,a,b,c)

    p = x^2+1
    @test int_polynomial_ellipsoid(p,a,b,c) ≈ 4/3*a*b*c + 4/3*a*b*c*(a^2*1/5)

    u = Mire.ptype{ComplexF64}[x^2,im*y,z]
    v = Mire.ptype{ComplexF64}[z^0,x,y]
    @test Mire.inner_product(u,v,cmat) == Mire.inner_product(u,v,a,b,c) ≈ 4/3*a*b*c*(a^2*1/5)

end


@testset "Projection functions" begin

    N = 7
    V = Ellipsoid(1.1,1.0,0.9)

    vbasis = LebovitzBasis(N,V).el
    nu = length(vbasis)
    A = spzeros(nu,nu)
    cmat = cacheint(N,V)

    # projectforce!(A,cmat,vbasis,vbasis,inertial)
    A = projectforce(vbasis,vbasis,cmat,inertial)

    # @test A ≈ B

    B = Mire.projectforcet(vbasis,vbasis,cmat,inertial)

    @test A ≈ B

    
end

@testset "Assemblies" begin
    Le = 1e-4
    Lu = 1e5
    Ω = [0,0,1.0]
    N = 5
    b0 = Mire.bpol_g21(1,1,0)
    prob = MHDProblem(N, Sphere(), Ω, Le, Lu, b0, QGIMBasis, InsulatingMFBasis)
    prob_t = deepcopy(prob)
    assemble!(prob)
    assemble!(prob_t; threads=true)
    
    @test prob.LHS ≈ prob_t.LHS
    @test prob.RHS ≈ prob_t.RHS

end