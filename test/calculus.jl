
@testset "Calculus function tests" begin


    p = x^2+y^2+2z^2
    @test Δ(p) == 8

    @test ∇(p) == [2x,2y,4z]

    u = [y*x,x^2*z*y,z^2]
    @test divergence(u) == ∇⋅u == y+x^2*z+2z

    @test curl(u) == ∇×u == [-x^2*y,0,2x*y*z-x]
    v = [0,0,x]

    @test advecterm(v,u) == (v⋅∇)(u) == [0,x^3*y,2x*z]

end
