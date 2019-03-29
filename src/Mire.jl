module Mire

using MultivariatePolynomials, TypedPolynomials, LinearAlgebra, SparseArrays, SpecialFunctions

export x,y,z,Π, ∇, Δ, div, curl, F, ex,ey,ez,
    combos, N1, N2, n_u, n_c,
    inertial, coriolis, lorentz, advection,
    eigen, vel, eigenvel, angularmom,r

include("assemble.jl")

export  assemblehd, assemblemhd, mat_force, mat_force_galerkin!, assemblemhd_diffusion

include("integration.jl")

export inner_product, int_monomial_ellipsoid, int_polynomial_ellipsoid, cacheint

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
advecterm(u,v) = [u[1]*∂(v[i],x) + u[2]*∂(v[i],y) + u[3]*∂(v[i],z) for i=1:3]
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
    @inbounds for k=0:N-1, i=0:N-1, j=0:N-1
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
lorentz(B,a,b,c,B0) = curl(B) × B0 + curl(B0) × B
# lorentz(B,a,b,c,B0) = curl(B) × B0 + curl(B0) × B + curl(B0) × B0 #experiment

advection(u,a,b,c,B0) = curl(u × B0)
diffusion(B,a,b,c,B0,Lu) = 1/Lu*Δ.(B+B0)

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
angularmom(u,cmat) = int_polynomial_ellipsoid((u×r)[3],cmat)

end #module
