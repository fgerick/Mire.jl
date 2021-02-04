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
    Le::T
    Lu::T
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
                    Le::T,
                    Lu::T,
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

    return MHDProblem(N, V, Ω, vptype{T}(B0), Le, Lu, vbasis, bbasis, cmat, LHS, RHS)
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
Le = 1/norm(Ω)

return MHDProblem(N, V, Ω, vptype{T}(B0), Le, T(Inf), vbasis, bbasis, cmat, LHS, RHS)
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
function assemble!(P::MHDProblem{T,V}; threads=false, kwargs...) where {T,V}
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu + nb

    # pfun! = threads ? projectforcet! : projectforce!
    if threads
        LHST = zeros(T,size(P.LHS)...)
        RHST = zeros(T,size(P.RHS)...)

        projectforcet!(view(LHST, 1:nu, 1:nu), P.cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t
        if typeof(P.bbasis) <: InsulatingMFBasis
            ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
            LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
            ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(ls)))
            projectforcet_symmetric_neighbours!(view(LHST,nu+1:nmat,nu+1:nmat),P.cmat,bbasis,bbasis, inertial,LS,MS,ispt; kwargs...) #∂j/∂t
        else
            projectforce!(view(LHST, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, inertial; kwargs...) #∂j/∂t
        end
        println("assemble LHS done!")

        projectforcet!(view(RHST, 1:nu, 1:nu), P.cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...) #Ω×u
        println("assemble Ω×u done!")
        projectforcet!(view(RHST, 1:nu, nu+1:nmat), P.cmat, vbasis, bbasis, lorentz, P.B0; kwargs...) #j×b
        println("assemble j×B done!")
        projectforcet!(view(RHST, nu+1:nmat, 1:nu), P.cmat, bbasis, vbasis, advection, P.B0; kwargs...)
        println("assemble ∇×u×B₀ done!")

        if !isinf(P.Lu)
            projectforcet!(view(RHST, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, b->1/P.Lu*diffusion(b); kwargs...)
            println("assemble 1/Lu ΔB² done!")
        end

        P.LHS = sparse(LHST)
        P.RHS = sparse(RHST)
    else
        projectforce!(view(P.LHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, inertial) #∂u/∂t
        projectforce!(view(P.LHS, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, inertial) #∂j/∂t
        projectforce!(view(P.RHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, coriolis, P.Ω) #Ω×u
        projectforce!(view(P.RHS, 1:nu, nu+1:nmat), P.cmat, vbasis, bbasis, lorentz, P.B0) #j×b
        projectforce!(view(P.RHS, nu+1:nmat, 1:nu), P.cmat, bbasis, vbasis, advection, P.B0)        
        if !isinf(P.Lu)
            projectforce!(view(P.RHS, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, b->1/P.Lu*diffusion(b); kwargs...)
        end
    end
    nothing
end


function normalizebasis!(P::MireProblem{T}; n_cache::Int = 10^6) where T
    ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache)
    el = map(u->u/sqrt(Mire.inner_product!(ptemp,u,u,P.cmat)) , P.vbasis.el)
    P.vbasis = typeof(P.vbasis)(P.N,P.V,el)

    if typeof(P) <: MHDProblem
       el = map(b->b/sqrt(Mire.inner_product!(ptemp,b,b,P.cmat)) , P.bbasis.el)
       P.bbasis = typeof(P.bbasis)(P.N,P.V,el)
    end
    return nothing
end
