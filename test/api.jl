
@testset "miscellaneous function tests" begin

    for n in 3:7
        vs = vel(n,1,1,1)
        @test length(vs) == n_u(n)
    end

    vs = [[x,y,z],ex,ey,ez]
    αs = [1//2,1,0,3]
    @test eigenvel(vs,αs) == [1//2*x+1,1//2*y,1//2*z+3]
    @test geo_vel(0,1,1,1)[1] == [-y,x,0]

end


@testset "monomial integrations" begin

    a,b,c = 0.2,0.3,0.4
    cmat = cacheint(3,a,b,c)

    @test cmat[1,1,1] ≈ int_monomial_ellipsoid(x^0,a,b,c) ≈ 4/3*a*b*c
    @test cmat[2,3,1] ≈ int_monomial_ellipsoid(x^2*y^3*z,a,b,c)

    p = x^2+1
    @test int_polynomial_ellipsoid(p,a,b,c) == 4/3*a*b*c + 4/3*a*b*c*(a^2*1/5)

    u = [x^2,im*y,z]
    v = [z^0,x,y]
    @test Mire.inner_product(cmat,u,v) == Mire.inner_product(u,v,a,b,c) ≈ 4/3*a*b*c*(a^2*1/5)
    @test Mire.inner_product_real(cmat,v,v) == Mire.inner_product(v,v,a,b,c) ≈ 4/3*a*b*c*(1+a^2*1/5+b^2*1/5)

end
