module Mire

using MultivariatePolynomials, TypedPolynomials, LinearAlgebra, SparseArrays, SpecialFunctions

export x,y,z,Π, ex,ey,ez,eigen, vel, eigenvel,
        assemblehd, assemblemhd,
        angularmom,r

@polyvar x y z


Π(n::Int,m::Int,l::Int) = x^n*y^m*z^l

Π(n::BigInt,m::BigInt,l::BigInt) = Π(Int(n),Int(m), Int(l))

∂ = differentiate
∇(ψ) = [∂.(ψ,(x,y,z))...]
Δ(ψ) = ∂(∂(ψ,x),x) + ∂(∂(ψ,y),y) + ∂(∂(ψ,z),z)
div(u) = ∂(u[1],x)+∂(u[2],y)+∂(u[3],z)
curl(u) = [∂(u[3],y)-∂(u[2],z),∂(u[1],z)-∂(u[3],x),∂(u[2],x)-∂(u[1],y)]

F(a,b,c) = (1-x^2/a^2-y^2/b^2-z^2/c^2)
const ex=[1,0,0]
const ey=[0,1,0]
const ez=[0,0,1]

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

function vel(N::Int,a,b,c)
    gp,hp = combos(N)
    v_1 = [v1(h...,a,b,c) for h in vcat(gp,hp)]
    v_2 = [v2(h...,a,b,c) for h in vcat(gp,hp)]
    v_3 = [v3(g...,a,b,c) for g in gp]

    return vcat(v_1,v_2,v_3)
end


N1(n) = n*(n+1)÷2
N2(n) = n*(n+1)*(n+2)÷6

n_c(N::Integer) = N1(n)+N2(N)
n_u(N::Int) = N1(N)+2N2(N)


## Force functions:

#hydro:
inertial(u,a,b,c) = u
coriolis(u,a,b,c,Ω) = -2*Ω×u
viscous(u,a,b,c,Ek) = Ek*Δ.(u)
viscous(u,a,b,c,Lu,Pm) = Pm/Lu*Δ.(u)

#magnetic:
lorentz(B,a,b,c,B0) = curl(B)×B0 + curl(B0)×B
advection(u,a,b,c,B0) = curl(u×B0)
diffusion(B,a,b,c,Lu) = 1/Lu*Δ.(B)

## create matrices using galerkin:

function mat_force_galerkin!(A::AbstractArray{T,2},vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A


    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...)
        for i=1:n_A
            A[i,j] = inner_product(vs[i],f,a,b,c)
        end
    end
end

function mat_force(N::Integer,vs, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin!(A,vs,N ,forcefun,a,b,c,args...)
    return A
end



### ellipsoid integration


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
        return 8big(π)*a^(2*γ₁+1)*b^(2*γ₂+1)*c^(2γ₃+1)*fg1*f2g/fg/f2g3
    else
        zero(BigFloat)
    end
end

int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))

inner_product(u,v,a::Real,b::Real,c::Real) = int_polynomial_ellipsoid(dot(u,v),a,b,c)



function assemblemhd(N,a,b,c,Ω,b0)
    T = typeof(a)
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    LHS = spzeros(T,2n_mat,2n_mat)
    RHS = spzeros(T,2n_mat,2n_mat)

    LHS[1:n_mat,1:n_mat] .= mat_force(N,vs,inertial,a,b,c)
    LHS[n_mat+1:end,n_mat+1:end] .= mat_force(N,vs,inertial,a,b,c)

    RHS[1:n_mat,1:n_mat] .= mat_force(N,vs,coriolis,a,b,c,Ω)
    RHS[1:n_mat,n_mat+1:end] .= mat_force(N,vs,lorentz,a,b,c,b0)

    RHS[n_mat+1:end,1:n_mat] .= mat_force(N,vs,advection,a,b,c,b0)

    return LHS,RHS, vs
end


function eigenvel(N::Integer,vs,λs,n_ev::Integer,a::T,b::T,c::T; norm =true) where {T<: Real, S<:Number}
    # ac=all_combos(N)
    # vs=[vel(ac[i]...,a,b,c) for i=1:length(ac)]
    # us = gramschmidt(vs,a,b,c)
    vo= sum([λs[i,n_ev]*vs[i] for i=1:length(vs)])
    return norm ? vo/√complex(inner_product(vo,vo,a,b,c)) : vo
end

const r = [x, y, z]
angularmom(u,a,b,c) = int_polynomial_ellipsoid((u×r)[3],a,b,c)

end #module
