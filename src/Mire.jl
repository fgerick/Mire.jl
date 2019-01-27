module Mire

using MultivariatePolynomials, TypedPolynomials, LinearAlgebra, SparseArrays, SpecialFunctions

export x,y,z,Π,poly_inertialforce,poly_coriolisforce,poly_lorentzforce, poly_magadvection,
        all_combos, mat_inertialforce, mat_coriolisforce,mat_inertialforce!,mat_coriolisforce!,
        mat_lorentzforce, mat_lorentzforce!, mat_maginertialforce, mat_maginertialforce!,
        mat_magadvection,mat_magadvection!,
        eigen, vel, eigenvel,
        assemblehd, assemblemhd

@polyvar x y z

# g(n::Int,m::Int) = x^n*y^m #*z^l
Π(n::Int,m::Int,l::Int) = x^n*y^m*z^l

# Π(n::BigInt,m::BigInt,l::BigInt) = Π(Int(n),Int(m), Int(l))

∂ = differentiate
∇(ψ) = [∂.(ψ,(x,y,z))...]
Δ(ψ) = ∂(∂(ψ,x),x) + ∂(∂(ψ,y),y) + ∂(∂(ψ,z),z)
div(u) = ∂(u[1],x)+∂(u[2],y)+∂(u[3],z)
curl(u) = [∂(u[3],y)-∂(u[2],z),∂(u[1],z)-∂(u[3],x),∂(u[2],x)-∂(u[1],y)]

F(a,b,c) = (1-x^2/a^2-y^2/b^2-z^2/c^2)
ex=[1,0,0]
ey=[0,1,0]
ez=[0,0,1]

v1(n::Int,m::Int,l::Int,a,b,c) = ∇(Π(n,m,l)*F(a,b,c))×ex
v2(n::Int,m::Int,l::Int,a,b,c) = ∇(Π(n,m,l)*F(a,b,c))×ey
v3(n::Int,m::Int,l::Int,a,b,c) = ∇(Π(n,m,l)*F(a,b,c))×ez



function combos(N::Int)
    gpairs = Vector{Int}[]
    hpairs = Vector{Int}[]
    @inbounds for i=0:N-1, j=0:N-1, k=0:N-1
        if (i + j + k <= N-1)
            if k==0
                # if i==j==0
                #     push!(hpairs,[i,j,k])
                # end
                push!(gpairs,[i,j,k])
            else
                push!(hpairs,[i,j,k])
            end
        end
    end
    return gpairs, hpairs
end

function v(N::Int,a,b,c)
    gp,hp = combos(N)
    v_1 = [v1(h...,a,b,c) for h in vcat(gp,hp)]
    v_2 = [v2(h...,a,b,c) for h in vcat(gp,hp)]
    v_3 = [v3(g...,a,b,c) for g in gp]

    return vcat(v_1,v_2,v_3)
end


N1(n) = n*(n+1)÷2
# nu(n)=n*(n+1)*(2n+7)÷6
N2(n) = n*(n+1)*(n+2)÷6
# n_u(n) = n*(n+1)÷2 + 2n*(n+1)*(n+2)÷6
n_c(N::Integer) = N1(n)+N2(N)
n_u(N::Int) = N1(N)+2N2(N)
# all_combos(N::Integer) = combos(N)


inertial(u,a,b,c) = u

coriolis(u,a,b,c,Ω) = -2*Ω×u

viscous(u,a,b,c,Ek) = Ek*Δ.(u)

function mat_force_galerkin!(A::AbstractArray{T,2},vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A


    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...)
        for i=1:n_A
            A[i,j] = inner_product(vs[i],f,a,b,c)
            # A[i,j] = inner_product_2D(p,Π(combos[j]...),a,b,c)
        end
    end
end

function mat_force(N::Integer,vs, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin!(A,vs,N ,forcefun,a,b,c,args...)
    # mat_force_coefficient!(A,N ,0,polyfun,a,b,c,args...)
    return A
end

function int_monomial_ellipsoid(p::Monomial,a::Real,b::Real,c::Real)
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    if iseven(i) && iseven(j) && iseven(k)
        γ₁ = i÷2
        γ₂ = j÷2
        γ₃ = k÷2
        γ = γ₁ + γ₂ + γ₃
        f2g = factorial(2γ₁)*factorial(2γ₂)*factorial(2γ₃)
        fg = factorial(γ₁)*factorial(γ₂)*factorial(γ₃)
        fg1 = factorial(big(γ₁+γ₂+γ₃+1))
        f2g3 = factorial(big(2γ₁+2γ₂+2γ₃+3))
        # f2g = factorial(2γ)
        # fg = factorial(γ₁)*factorial(γ₂)*factorial(γ₃)
        # fg1 = factorial(γ+1)
        # f2g3 = factorial(2γ+3)
        return 8big(π)*a^(2*γ₁+1)*b^(2*γ₂+1)*c^(2γ₃+1)*fg1*f2g/fg/f2g3
    else
        zero(BigFloat)
    end
end

int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))

inner_product(u,v,a::Real,b::Real,c::Real) = int_polynomial_ellipsoid(dot(u,v),a,b,c)



end #module
