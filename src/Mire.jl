module Mire

using MultivariatePolynomials, TypedPolynomials, LinearAlgebra, SparseArrays, SpecialFunctions

export x, y, z, r, Π, ∇, Δ, gradient, laplacian, divergence, curl, advecterm
export F, ex,ey,ez, combos, N1, N2, n_u, n_c
export inertial, coriolis, lorentz, advection
export vel, eigenvel, qg_vel, geo_vel, velocities, magneticfields

include("assemble.jl")

export  assemblehd, assemblemhd, assemblemhd!, projectforce, projectforce!

include("integration.jl")

export inner_product, int_monomial_ellipsoid, int_polynomial_ellipsoid, cacheint, cacheint_surface_torque

# Cartesian coordinates as polynomial variables

@polyvar x y z

# Monomials
Π(n::Integer,m::Integer,l::Integer) = x^n*y^m*z^l

# Calculus definitions
∂ = differentiate
gradient(ψ) = [∂.(ψ,(x,y,z))...]
∇ = gradient
laplacian(ψ) = ∂(∂(ψ,x),x) + ∂(∂(ψ,y),y) + ∂(∂(ψ,z),z)
Δ = laplacian
divergence(u) = ∂(u[1],x)+∂(u[2],y)+∂(u[3],z)
curl(u) = [∂(u[3],y)-∂(u[2],z),∂(u[1],z)-∂(u[3],x),∂(u[2],x)-∂(u[1],y)]
advecterm(u,v) = [u[1]*∂(v[i],x) + u[2]*∂(v[i],y) + u[3]*∂(v[i],z) for i=1:3]

# Lebovitz 1989, eq. (39b)
F(a,b,c) = (1-x^2/a^2-y^2/b^2-z^2/c^2)

# Cartesian unit vectors
const ex = [1,0,0]
const ey = [0,1,0]
const ez = [0,0,1]

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

# Number of basis vectors (Lebovitz 1989/Vidal, Cebron 2017).
N1(n::Int) = n*(n+1)÷2
N2(n::Int) = n*(n+1)*(n+2)÷6

n_u(N::Int) = N1(N)+2N2(N)


## Force functions:

#hydro:
inertial(u) = u
coriolis(u,Ω)  = -2*Ω×u

#magnetic:
lorentz(B,B0)  = curl(B) × B0 + curl(B0) × B
advection(u,B0)  = curl(u × B0)

## reconstruction of velocity from eigenvector

"""
    eigenvel(v,α)
Reconstructs velocity u=∑αᵢvᵢ
"""
function eigenvel(v,α)
    @assert length(v)==length(α) "Coefficients should have the same length as basis vectors"
    sum([α[i]*v[i] for i=1:length(v)])
end

eigenvel(vs,αs,n_ev::Integer) = eigenvel(vs,αs[:,n_ev])

function velocities(vs,αs)
    nv = length(vs)
    [eigenvel(vs,αs[1:nv,i]) for i=1:size(αs,2)]
end

function magneticfields(bs,αs)
    nb = length(bs)
    [eigenvel(bs,αs[end-nb+1:end,i]) for i=1:size(αs,2)]
end


## QG tools
qg_combos(N::Integer) = [[i,j] for i=0:N for j=0:N if (i+j<N)]
# qg_combos(N::Integer) = (N==0) ? [[0,0]] : vcat(combos(N-1),[[i,j] for i=0:N for j=0:N if (N-1<i+j<=N)]) #sorted

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

"""
    qg_vel(N::Integer,a::T,b::T,c::T) where T

Generate QG basis vectors.
"""
function qg_vel(N::Int,a::T,b::T,c::T) where T
    cs = qg_combos(N)
    return [qg_vel(ci...,a,b,c) for ci in cs]
end

"""
    geo_veln(n::Integer,a::T,b::T,c::T) where T

Generate geostrophic basis vector \$\\mathbf{u}_n=\\nabla(h^{3+2n})\\times\\mathbf{1}_z\$
"""
function geo_veln(n::Integer,a::T,b::T,c::T) where T
    h2 = c^2*(1-x^2/a^2-y^2/b^2)
    ez = [0,0,1]
    hgradh = [-c^2*x/a^2,-c^2*y/b^2,0]
    return (3+2n)//3*h2^n*z^0 * hgradh×ez
end

"""
    geo_vel(N::Integer,a::T,b::T,c::T) where T

Generate geostrophic basis vectors.
"""
function geo_vel(N::Int,a::T,b::T,c::T) where T
    return [geo_veln(n,a,b,c) for n in 0:N]
end


#Quagmire polynomial functions (2D reduced set of equations)

Π(n::Integer, m::Integer) = x^n*y^m
"""
    poly_inertial(n::Integer,m::Integer,a::T,b::T,c::T) where T

Polynomial for basis (n,m) of inertial force
"""
function poly_inertial(n::Integer,m::Integer,a::T,b::T,c::T) where T
    out = -((n+3)*(n+1)*c^2/a^2+(m+3)*(m+1)*c^2/b^2)*Π(n,m)+c^2*n*(n-1)*Π(n-2,m)
    out += c^2*m*(m-1)*Π(n,m-2)-n*(n-1)*c^2/b^2*Π(n-2,m+2)-m*(m-1)*c^2/a^2*Π(n+2,m-2)
    out += c^4/3*(n*(n-1)/b^4*Π(n-2,m+2) + m*(m-1)/a^4*Π(n+2,m-2)-(2n*m+n+m)/(a^2*b^2)*Π(n,m))
    return out
end

"""
    poly_coriolis(n::Integer,m::Integer,a::T,b::T,c::T,Ω::T) where T

Polynomial for basis (n,m) of Coriolis force.
"""
function poly_coriolis(n::Integer,m::Integer,a::T,b::T,c::T,Ω::T) where T
    E = -c^2/2*(x^2/a^2+y^2/b^2-1)
    ez = [0,0,1]
    return -2Ω*(∇(E) × ∇(Π(n,m)))⋅ez
end

#mhd:
"""
    poly_lorentz(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T

Polynomial for basis (n,m) of Lorentz force.
"""
function poly_lorentz(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T
    coeffs_A0, nm_pairs_A0 = coefficients(A0), exponents.(monomials(A0))
    D2A0 = sum([cf*poly_inertial(nmpair...,a,b,c) for (cf,nmpair) in zip(coeffs_A0,nm_pairs_A0)])
    D2A = poly_inertial(n,m,a,b,c)
    h2 = -c^2*(x^2/a^2+y^2/b^2-1)
    E = h2/2
    ez = [0,0,1]
    out =  h2.*(∇(D2A0)×∇(Π(n,m)))⋅ez .+ 3*Π(n,m) .*(∇(D2A0)×∇(E))⋅ez  .- D2A0.*(∇(E)×∇(Π(n,m)))⋅ez
    out += h2.*(∇(D2A) ×∇(A0))⋅ez     .+ 3*A0     .*(∇(D2A) ×∇(E))⋅ez  .- D2A .*(∇(E)×∇(A0))⋅ez
    return out
end

"""
    poly_maginertial(n::Integer,m::Integer,a::T,b::T,c::T) where T


"""
function poly_inertialmag(n::Integer,m::Integer,a::T,b::T,c::T) where T
    return Π(n,m)
end

"""
    poly_magadvection(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T


"""
function poly_advection(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T
    h2 = -c^2*(x^2/a^2+y^2/b^2-1)
    E = h2/2
    ez = [0,0,1]
    out = h2.*(∇(Π(n,m)) × ∇(A0)) ⋅ ez .+ 3 .*A0 .*(∇(Π(n,m))×∇(E))⋅ez .+ 3 .*Π(n,m) .*(∇(E)×∇(A0))⋅ez
    return out
end





end #module
