### Functions to create matrices and assemble the full system.

abstract type MireProblem{T,V} end


"""
    HDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T,Vol}

Defines hydrodynamic problem.

Example:

```
N = 5
Ω = [0,0,1.0]
V = Ellipsoid(1.1,1.0,0.9)
problem = HDProblem(N,V,Ω,LebovitzBasis)
```
"""
mutable struct HDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T,Vol}
    N::Int
    V::Vol
    Ω::Vector{T}
    vbasis::VectorBasis{T,Vol}
    cmat::Array{T,3}
    LHS::AbstractMatrix{T}
    RHS::AbstractMatrix{T}
end

"""
    MHDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T,Vol}

Defines magnetohydrodynamic problem.

Example for a hybrid QG model with 3-D magnetic field with 
conducting boundary condition and QG velocities:

```
N = 5
Ω = [0,0,1.0]
a,b,c = 1.1,1.0,0.9
V = Ellipsoid(a,b,c)
B0 = [-y/b^2,x/a^2,0] #Malkus field
problem = MHDProblem(N,V,Ω,B0,QGBasis,ConductingMFBasis)
```
"""
mutable struct MHDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T,Vol}
    N::Int
    V::Vol
    Ω::Vector{T}
    B0::vptype{T}
    vbasis::VectorBasis{T,Vol}
    bbasis::VectorBasis{T,Vol}
    cmat::Array{T,3}
    LHS::AbstractMatrix{T}
    RHS::AbstractMatrix{T}
end

function HDProblem(
                    N::Int,
                    V::Volume{T},
                    Ω::Vector{T},
                    ::Type{VB},
                ) where {T<:Number,VB<:VectorBasis}

    vbasis = VB(N, V)
    cmat = cacheint(N, V)
    n = length(vbasis.el)
    LHS = spzeros(T, n, n)
    RHS = spzeros(T, n, n)

    return HDProblem(N, V, Ω, vbasis, cmat, LHS, RHS)
end

function MHDProblem(
                    N::Int,
                    V::Volume{T},
                    Ω::Vector{T},
                    B0,
                    ::Type{VB},
                    ::Type{BB},
                ) where {T<:Number,VB<:VectorBasis,BB<:VectorBasis}

    vbasis = VB(N, V)
    bbasis = BB(N, V)
    cmat = cacheint(N, V)
    nu = length(vbasis.el)
    nb = length(bbasis.el)
    n = nu + nb
    LHS = spzeros(T, n, n)
    RHS = spzeros(T, n, n)

    return MHDProblem(N, V, Ω, vptype{T}(B0), vbasis, bbasis, cmat, LHS, RHS)
end

"""
    assemble!(P::HDProblem{T,V}) where {T,V}

Assembles the matrices `P.LHS` and `P.RHS`, i.e. projecting the velocity basis
`P.vbasis` on the inertial acceleration and Coriolis force.
"""
function assemble!(P::HDProblem{T,V}) where {T,V}

    projectforce!(P.LHS, P.cmat, P.vbasis.el, inertial)
    projectforce!(P.RHS, P.cmat, P.vbasis.el, coriolis, P.Ω)

    return nothing
end

"""
    assemble!(P::MHDProblem{T,V}) where {T,V}

Assembles the matrices `P.LHS` and `P.RHS`, i.e. projecting the velocity and
magnetic field bases on the inertial acceleration, Coriolis force, Lorentz force
and mgnetic advection.
"""
function assemble!(P::MHDProblem{T,V}) where {T,V}
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu + nb
    projectforce!(view(P.LHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, inertial) #∂u/∂t
    projectforce!(view(P.LHS, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, inertial) #∂j/∂t
    projectforce!(view(P.RHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, coriolis, P.Ω) #Ω×u
    projectforce!(view(P.RHS, 1:nu, nu+1:nmat), P.cmat, vbasis, bbasis, lorentz, P.B0) #j×b
    projectforce!(view(P.RHS, nu+1:nmat, 1:nu), P.cmat, bbasis, vbasis, advection, P.B0)
    nothing
end