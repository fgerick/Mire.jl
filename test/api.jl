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
    @test int_polynomial_ellipsoid(p,a,b,c) == 4/3*a*b*c + 4/3*a*b*c*(a^2*1/5)

    u = Mire.ptype{ComplexF64}[x^2,im*y,z]
    v = Mire.ptype{ComplexF64}[z^0,x,y]
    @test Mire.inner_product(cmat,u,v) == Mire.inner_product(u,v,a,b,c) ≈ 4/3*a*b*c*(a^2*1/5)
    @test Mire.inner_product_real(cmat,v,v) == Mire.inner_product(v,v,a,b,c) ≈ 4/3*a*b*c*(1+a^2*1/5+b^2*1/5)

end

@testset "Insulating magnetic field basis" begin
    N = 3
    V = Ellipsoid(3//5,4//5,6//5)
    # check if ellipsoid throws an error:
    @test_throws ArgumentError InsulatingMFBasis(N,V; norm=Mire.Nonorm{Rational{Int}}())

    #first element of basis check:
    b = InsulatingMFBasis(N, Sphere{Rational{Int}}(), norm=Mire.Nonorm{Rational{Int}}())
    l, m, n = 1, -1, 0
    RLM = - y
    P1m1 = -(x^2+y^2+z^2)^(n + 1)*RLM*1//(2(n+1)*(2l+2n+3))
    ∇vi = ∇(-Mire.vi(l,n)*RLM)
    
    @test b.el[1] == ∇ × (∇ × (P1m1*[x,y,z])) + ∇vi

    # Numerically check if any current normal to the boundary exists:
    j = curl.(b.el)

    function normalmax(u)
        m = 0.0
        dotn=dot(u,[x,y,z])
        @inbounds @fastmath for theta in 0:0.1:pi, phi in 0:0.1:2pi
            m = max(m,abs(dotn(x=>cos(phi)*sin(theta),y=>sin(phi)*sin(theta),z=>cos(theta))))
        end
        return m
    end
    
    jn_at_surface = normalmax.(j)

    @test all( isapprox.(jn_at_surface,0, atol=1e-13))

    #test LMN function
    ls,ms,ns,lstor,mstor,nstor = Mire.LMN(b)
    @test maximum(ls)==N
    @test ms[1] == -1
    @test mstor[1] == -1
    @test maximum(mstor) == N-1

    ## todo: Bpol + ∇Φi = ∇Φe at r = 1.

end

@testset "Projection functions" begin

    N = 7
    V = Ellipsoid(1.1,1.0,0.9)

    vbasis = LebovitzBasis(N,V).el
    nu = length(vbasis)
    A = spzeros(nu,nu)
    cmat = cacheint(N,V)

    projectforce!(A,cmat,vbasis,vbasis,inertial)
    B = projectforce(cmat,vbasis,vbasis,inertial)

    @test A ≈ B

    C = zeros(nu,nu)
    Mire.projectforcet!(C,cmat,vbasis,vbasis,inertial)

    @test B ≈ sparse(C)

    
end