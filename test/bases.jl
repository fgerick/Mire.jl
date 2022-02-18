
function normalmax(u,a,b,c)
    m = 0.0
    dotn=dot(u,[x/a^2,y/b^2,z/c^2])
    @inbounds @fastmath for theta in 0:0.1:pi, phi in 0:0.1:2pi
        m = max(m,abs(dotn(x=>a*cos(phi)*sin(theta),y=>b*sin(phi)*sin(theta),z=>c*cos(theta))))
    end
    return m
end
    
@testset "QGBasis" begin

    a,b,c = 1,2,3
    basis = QGBasis(7,Ellipsoid(a,b,c))
    @test ∇⋅(sum(basis.el)) ≈ 0
    @test normalmax(sum(basis.el),a,b,c) ≈ 0 atol=sqrt(eps()) 

end

@testset "QGIMBasis" begin

    basis = QGIMBasis(7,Sphere())
    cmat = cacheint(7,Sphere())*pi
    @test ∇⋅(sum(basis.el)) ≈ 0
    @test normalmax(sum(basis.el),ones(3)...) ≈ 0 atol=sqrt(eps()) 

    #check orthogonality:
    @test [inner_product(ui,uj,cmat) for ui in basis.el, uj in basis.el] ≈ I 
end

@testset "QGRIMBasis" begin

    basis = QGRIMBasis(7,Sphere())
    @test ∇⋅(sum(basis.el)) ≈ 0
    @test normalmax(sum(basis.el),ones(3)...) ≈ 0 atol=sqrt(eps()) 
end


@testset "LebovitzBasis" begin
    a,b,c = 1,2,3
    @test Mire.v1(0,0,0,Ellipsoid(a,b,c)) == [0, -2//1*z/c^2, 2//1*y/b^2]
    @test Mire.v2(0,0,0,Ellipsoid(a,b,c)) == [2//1*z/c^2, 0, -2//1*x/a^2]
    @test Mire.v3(0,0,0,Ellipsoid(a,b,c)) == [-2//1*y/b^2, 2//1*x/a^2, 0]
     
    basis = LebovitzBasis(7,Ellipsoid(a,b,c))
    @test ∇⋅(sum(basis.el)) ≈ 0
    @test normalmax(sum(basis.el),a,b,c) ≈ 0 atol=sqrt(eps())
end


@testset "Insulating magnetic field basis" begin
    N = 3
    V = Ellipsoid(3//5,4//5,6//5)
    # check if ellipsoid throws an error:
    @test_throws MethodError InsulatingMFBasis(N,V; norm=Mire.Nonorm{Rational{Int}})

    #first element of basis check:
    b = InsulatingMFBasis(N, Sphere{Rational{Int}}(), norm=Mire.Nonorm{Rational{Int}})
    l, m, n = 1, -1, 0
    RLM = - y
    P1m1 = -(x^2+y^2+z^2)^(n + 1)*RLM*1//(2(n+1)*(2l+2n+3))
    ∇vi = ∇(-Mire.vi(l,n)*RLM)
    
    @test b.el[1] == ∇ × (∇ × (P1m1*[x,y,z])) + ∇vi

    # Numerically check if any current normal to the boundary exists:
    j = curl.(b.el)

    
    jn_at_surface = normalmax.(j,ones(3)...)

    @test all( isapprox.(jn_at_surface,0, atol=1e-13))

   
    @test maximum(map(u->maximum(maxdegree.(u)),b.el)) == N 
    #test LMN function
    # ls,ms,ns,lstor,mstor,nstor = Mire.LMN(b)
    # @test maximum(lstor)==N-1
    # @test mstor[1] == -1
    # @test ms[1] == -1
    # @test maximum(ms) == N-2

    ## todo: Bpol + ∇Φi = ∇Φe at r = 1.

end

@testset "Ivers basis" begin
    
    Ω = [0,0,1]
    V = Ellipsoid(1.3,0.9,1.1)
    N = 3
    evals = []

    prob1 = HDProblem(N,V, Ω, LebovitzBasis)
    assemble!(prob1)
    esol1 = eigen(Matrix(prob1.RHS), Matrix(prob1.LHS))
    ω1 = esol1.values

    prob2 = HDProblem(N,V, Ω, IversBasis)
    assemble!(prob2)
    esol2 = eigen(Matrix(prob2.RHS), Matrix(prob2.LHS))
    ω2 = esol2.values

    @test maximum(map(u->maximum(maxdegree.(u)),prob2.vbasis.el)) == N
    @test sort(abs.(ω1)) ≈ sort(abs.(ω2))

end