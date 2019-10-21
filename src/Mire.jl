module Mire

using MultivariatePolynomials, TypedPolynomials, LinearAlgebra, SparseArrays, SpecialFunctions, HCubature

export x,y,z,Π, ∇, Δ, div, curl, F, ex,ey,ez,
    combos, N1, N2, n_u, n_c,
    inertial, coriolis, lorentz, advection,
    eigen, vel, eigenvel,r

include("assemble.jl")

export  assemblehd, assemblemhd, assemblemhd_hybrid, projectforce, projectforce!

include("integration.jl")

export inner_product, int_monomial_ellipsoid, int_polynomial_ellipsoid, cacheint, cacheint_surface_torque


# include miscellaneous function for analysis
# include("misc/analysis.jl")
# include("misc/tracking.jl")
#don't export them -> need Mire. prefix

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
    [[i,j,k] for i=0:N for j=0:N for k=0 if (i+j+k<N)],
    [[i,j,k] for i=0:N for j=0:N for k=1:N if (i+j+k<N)]
 end


"""
    vel(N,a,b,c)

Compute all velocity basis vectors for a given maximal degree N (Lebovitz 1989, eq. 41-42)
"""
function vel(N::Int,a::T,b::T,c::T) where T
    gp,hp = combos(N)
    v_1 = [v1(h...,a,b,c) for h in vcat(gp,hp)]
    v_2 = [v2(h...,a,b,c) for h in vcat(gp,hp)]
    v_3 = [v3(g...,a,b,c) for g in gp]

    return vcat(v_1,v_2,v_3)
end

# Number of basis vectors (Lebovitz/Vidal, Cebron 2017).
N1(n::Int) = n*(n+1)÷2
N2(n::Int) = n*(n+1)*(n+2)÷6

n_c(N::Int) = N1(n)+N2(N)
n_u(N::Int) = N1(N)+2N2(N)


## Force functions:

#hydro:
inertial(u::Array{P,1}) where {T, P<:Polynomial{T}} = u
coriolis(u::Array{P,1},Ω) where {T, P<:Polynomial{T}}  = -2*Ω×u

#magnetic:
lorentz(B::Array{P,1},B0) where {T, P<:Polynomial{T}}  = curl(B) × B0 + curl(B0) × B

inertialmag(B::Array{P,1}) where {T, P<:Polynomial{T}}  = B
advection(u::Array{P,1},B0) where {T, P<:Polynomial{T}}  = curl(u × B0)

"""
    eigenvel(v,α)
Reconstructs velocity u=∑αᵢvᵢ
"""
function eigenvel(v,α)
    @assert length(v)==length(α) "Coefficients should have the same length as basis vectors"
    sum([α[i]*v[i] for i=1:length(v)])
end

eigenvel(vs,αs,n_ev::Integer) = eigenvel(N,vs,αs[:,n_ev])


#QG tools
qg_combos(N::Integer) = [[i,j] for i=0:N for j=0:N if (i+j<=N)]

"""
    qg_vel(n::Integer,m::Integer,a::T,b::T,c::T) where T

Generate basis vector \$\\mathbf{u}=\\nabla(h^3x^ny^m)\\times\\nabla(z/h)\$
"""
function qg_vel(n::Integer,m::Integer,a::T,b::T,c::T) where T
    h2 = c^2*(1-x^2/a^2-y^2/b^2)
    ez = [0,0,1]
    hgradh = [-c^2*x/a^2,-c^2*y/b^2,0]
    return h2*∇(x^n*y^m)×ez+3*x^n*y^m*hgradh×ez-z*∇(x^n*y^m)×hgradh
end

function qg_vel(N::Int,a::T,b::T,c::T) where T
    cs=qg_combos(N)
    return [qg_vel(ci...,a,b,c) for ci in cs]
end


end #module
