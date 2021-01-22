module Mire

using CartesianSphericalHarmonics
using LinearAlgebra
using MultivariatePolynomials
using SparseArrays
using SpecialFunctions
using TypedPolynomials

export x, y, z, r, ∇, Δ, gradient, laplacian, divergence, curl, advecterm
export inertial, coriolis, lorentz, advection
export vel, eigenvel, velocities, magneticfields

# Cartesian coordinates as polynomial variables
const x = TypedPolynomials.Variable{:x}()
const y = TypedPolynomials.Variable{:y}()
const z = TypedPolynomials.Variable{:z}()

# different bases functions

include("bases.jl")
export Volume, Ellipsoid, Sphere, LebovitzBasis, QGBasis, ConductingMFBasis

include("assemble.jl")
export MireProblem, HDProblem, MHDProblem

include("projectfuns.jl")
export assemble!, projectforce, projectforce!

include("integration.jl")
export inner_product, int_monomial_ellipsoid, int_polynomial_ellipsoid, cacheint, cacheint_surface_torque


# Calculus definitions
∂ = differentiate
gradient(ψ) = [∂.(ψ,(x,y,z))...]
∇ = gradient
laplacian(ψ) = ∂(∂(ψ,x),x) + ∂(∂(ψ,y),y) + ∂(∂(ψ,z),z)
Δ = laplacian
divergence(u) = ∂(u[1],x)+∂(u[2],y)+∂(u[3],z)
curl(u) = [∂(u[3],y)-∂(u[2],z),∂(u[1],z)-∂(u[3],x),∂(u[2],x)-∂(u[1],y)]
advecterm(u,v) = [u[1]*∂(v[i],x) + u[2]*∂(v[i],y) + u[3]*∂(v[i],z) for i=1:3]

# Cartesian unit vectors
const ex = [1,0,0]
const ey = [0,1,0]
const ez = [0,0,1]

const r = [x, y, z]

## Force functions:

#hydro:
inertial(u) = u
coriolis(u,Ω)  = -2*Ω×u

#magnetic:
lorentz(B,B0)  = curl(B) × B0 + curl(B0) × B
advection(u,B0)  = curl(u × B0)
diffusion(B) = Δ.(B)



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



#2D reduced set, submodule
include("quagmire.jl")

end #module
