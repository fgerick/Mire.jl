module Mire

using MultivariatePolynomials, TypedPolynomials, LinearAlgebra, SparseArrays, SpecialFunctions

export x,y,z,Π, ex,ey,ez,eigen, vel, eigenvel,
        assemblehd, assemblemhd,
        angularmom,r

# Cartesian coordinates as polynomial variables
@polyvar x y z

# Monomials
Π(n::Int,m::Int,l::Int) = x^n*y^m*z^l
Π(n::BigInt,m::BigInt,l::BigInt) = Π(Int(n),Int(m), Int(l))

# Calculus definitions
∂ = differentiate
∇(ψ) = [∂.(ψ,(x,y,z))...]
Δ(ψ) = ∂(∂(ψ,x),x) + ∂(∂(ψ,y),y) + ∂(∂(ψ,z),z)
div(u) = ∂(u[1],x)+∂(u[2],y)+∂(u[3],z)
curl(u) = [∂(u[3],y)-∂(u[2],z),∂(u[1],z)-∂(u[3],x),∂(u[2],x)-∂(u[1],y)]

# Lebovitz 1989, eq. (39b)
F(a,b,c) = (1-x^2/a^2-y^2/b^2-z^2/c^2)

# Cartesian unit vectors
const ex=[1,0,0]
const ey=[0,1,0]
const ez=[0,0,1]

const r = [x, y, z]


# Basis vectors, Lebovitz 1989, eq. (41-42)
v1(n::Int,m::Int,l::Int,a,b,c) = ∇(Π(n,m,l)*F(a,b,c))×ex
v2(n::Int,m::Int,l::Int,a,b,c) = ∇(Π(n,m,l)*F(a,b,c))×ey
v3(n::Int,m::Int,l::Int,a,b,c) = ∇(Π(n,m,l)*F(a,b,c))×ez


"""
    combos(N::Int)

Computes all combos i,j,k satisfying i+j+k ≤ N.
"""
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


"""
    vel(N,a,b,c)

Compute all velocity basis vectors for a given maximal degree N (Lebovitz 1989, eq. 41-42)
"""
function vel(N::Int,a,b,c)
    gp,hp = combos(N)
    v_1 = [v1(h...,a,b,c) for h in vcat(gp,hp)]
    v_2 = [v2(h...,a,b,c) for h in vcat(gp,hp)]
    v_3 = [v3(g...,a,b,c) for g in gp]

    return vcat(v_1,v_2,v_3)
end

# Number of basis vectors (Lebovitz/Vidal, Cebron 2017).
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

"""
    mat_force_galerkin!(A,vs,N,forcefun,a,b,c, args...)

Fills Matrix `A` with Galerkin coefficients of force given by the function `forcefun(u,a,b,c,args...)`.
"""
function mat_force_galerkin!(A::AbstractArray{T,2},vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A


    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            A[i,j] = inner_product(vs[i],f,a,b,c) # calculates ∫ <uᵢ,f(uⱼ)> dV
        end
    end
end


function mat_force_galerkin_cached!(A::AbstractArray{T,2},cmat,vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A


    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            # A[i,j] = inner_product(vs[i],f,a,b,c) # calculates ∫ <uᵢ,f(uⱼ)> dV
            A[i,j] = inner_product_cached(cmat,vs[i],f)
        end
    end
end

"""
    mat_force(N,vs,forcefun,a,b,c, args...)

Allocates new matrix `A` and fills elements by calling
mat_force_galerkin!(A,vs,N,forcefun,a,b,c, args...).
"""
function mat_force(N::Integer,vs, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin!(A,vs,N ,forcefun,a,b,c,args...)
    return A
end

function mat_force_cached(N::Integer,cmat,vs, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin_cached!(A,cmat,vs,N ,forcefun,a,b,c,args...)
    return A
end




### ellipsoid integration

"""
Integral over the surface of an ellipsoid.
"""
function int_ellipsoid_surface(p::Monomial,a::Real,b::Real,c::Real)
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    if iseven(i) && iseven(j) && iseven(k)

        f1 = factorial(i)*factorial(j)*factorial(k)
        f2 = factorial(i÷2)*factorial(j÷2)*factorial(k÷2)
        f3 = factorial(2+i+j+k)
        f4= factorial(1+(i+j+k)÷2)
        ft = f1/f2/f3*f4
        return √big(π)*a^(i+1)*b^(j+1)*c^(k+1)*ft
    else
        zero(BigFloat)
    end
end

"""
    int_monomial_ellipsoid(p,a,b,c)

Integrate monomial `p=xⁱyʲzᵏ` over ellipsoid of semi-axes `a,b,c`.
"""
function int_monomial_ellipsoid(p::Monomial,a::Real,b::Real,c::Real)
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    return int_monomial_ellipsoid(i,j,k,a,b,c)
end

"""
    int_monomial_ellipsoid(i,j,k,a,b,c)

Integrate monomial `xⁱyʲzᵏ` over ellipsoid of semi-axes `a,b,c`.
"""
function int_monomial_ellipsoid(i::BigInt,j::BigInt,k::BigInt,a::Real,b::Real,c::Real)
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

int_polynomial_ellipsoid_surface(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_ellipsoid_surface.(monomial.(terms(p)),a,b,c))

"""
    inner_product(u,v,a,b,c)

Defines inner product in an ellipsoidal volume \$\\int\\langle u,v\\rangle dV\$.
"""
inner_product(u,v,a::Real,b::Real,c::Real) = int_polynomial_ellipsoid(dot(u,v),a,b,c)

function inner_product_cached(cmat,u,v)
    duv = dot(u,v)
    ip = zero(eltype(cmat))
    cs = coefficients(duv)
    exps = exponents.(monomial.(terms(duv)))
    @inbounds @simd for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end
"""
Function to precalculate monomial integrations.
"""
function cacheint(n::Int,a::T,b::T,c::T) where T<:Real
    Nmax=4n
    cachedmat=zeros(T,Nmax+1,Nmax+1,Nmax+1)
    for i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid(big(i),big(j),big(k),a,b,c)
    end
    return cachedmat
end


"""
    assemblemhd(N,a,b,c,Ω,b0)

Assembles MHD eigen system, such that
λAx=Bx

This is the dissipationless model, with

∂ₜu = -2Ω×u + (∇×b0)×b + (∇×b)×b0
∂ₜb = ∇×(u×b0)

with Ω = 1/Le * eΩ.
"""
function assemblemhd(N,a,b,c,Ω,b0)
    T = typeof(a)
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    A[1:n_mat,1:n_mat] .= mat_force(N,vs,inertial,a,b,c)
    A[n_mat+1:end,n_mat+1:end] .= mat_force(N,vs,inertial,a,b,c)

    B[1:n_mat,1:n_mat] .= mat_force(N,vs,coriolis,a,b,c,Ω)
    B[1:n_mat,n_mat+1:end] .= mat_force(N,vs,lorentz,a,b,c,b0)

    B[n_mat+1:end,1:n_mat] .= mat_force(N,vs,advection,a,b,c,b0)

    return A,B, vs
end
function assemblemhd_cachedint(N,cmat,a,b,c,Ω,b0)
    T = typeof(a)
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    A[1:n_mat,1:n_mat] .= mat_force_cached(N,cmat,vs,inertial,a,b,c)
    A[n_mat+1:end,n_mat+1:end] .= mat_force_cached(N,cmat,vs,inertial,a,b,c)

    B[1:n_mat,1:n_mat] .= mat_force_cached(N,cmat,vs,coriolis,a,b,c,Ω)
    B[1:n_mat,n_mat+1:end] .= mat_force_cached(N,cmat,vs,lorentz,a,b,c,b0)

    B[n_mat+1:end,1:n_mat] .= mat_force_cached(N,cmat,vs,advection,a,b,c,b0)

    return A,B, vs
end


"""
    eigenvel(N,vs,αs,a,b,c; norm=true)
Reconstructs velocity u, following Vidal & Cebron 2017 eq. (3.5).

If `norm` keyword is set `true` the velocity is normalised to satisfy \$\\int\\langle u,v\\rangle dV\$.
"""
function eigenvel(N::Integer,vs,αs,a::T,b::T,c::T; norm =true) where T<: Real
    vo= sum([αs[i]*vs[i] for i=1:length(vs)])
    return norm ? vo/√complex(inner_product(vo,vo,a,b,c)) : vo
end

eigenvel(N::Integer,vs,αs,n_ev::Integer,a::T,b::T,c::T; norm =true) where T<: Real = eigenvel(N,vs,αs[:,n_ev],a,b,c,norm=norm)

"""
    angularmom(u,a,b,c)

Calculates z-component of angular momentum.
"""
angularmom(u,a,b,c) = int_polynomial_ellipsoid((u×r)[3],a,b,c)

end #module
