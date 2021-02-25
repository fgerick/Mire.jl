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
    LHS::Union{AbstractMatrix{T},AbstractMatrix{Complex{T}}}
    RHS::Union{AbstractMatrix{T},AbstractMatrix{Complex{T}}}
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
    Le::T
    Lu::T
    B0::Union{vptype{T},vptype{Complex{T}}}
    vbasis::VectorBasis{T,Vol}
    bbasis::VectorBasis{T,Vol}
    cmat::Array{T,3}
    LHS::Union{AbstractMatrix{T},AbstractMatrix{Complex{T}}}
    RHS::Union{AbstractMatrix{T},AbstractMatrix{Complex{T}}}
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
    TM = coefficienttype(vbasis.el[1][1])
    LHS = spzeros(TM, n, n)
    RHS = spzeros(TM, n, n)

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
    TM = promote_type(coefficienttype(vbasis.el[1][1]),coefficienttype(bbasis.el[1][1]))
    LHS = spzeros(TM, n, n)
    RHS = spzeros(TM, n, n)
    Ω /= Le*norm(Ω)

    return MHDProblem(N, V, Ω, Le, Lu, vptype{TM}(B0), vbasis, bbasis, cmat, LHS, RHS)
end

function MHDProblem(
    N::Int,
    V::Volume{T},
    Ω::Vector{T},
    Le::T,
    B0,
    ::Type{VB},
    ::Type{BB},
) where {T<:Number,VB<:VectorBasis,BB<:VectorBasis}

    return MHDProblem(N, V, Ω, Le, T(Inf), vptype{T}(B0), VB, BB)
end

"""
    assemble!(P::HDProblem{T,V}) where {T,V}

Assembles the matrices `P.LHS` and `P.RHS`, i.e. projecting the velocity basis
`P.vbasis` on the inertial acceleration and Coriolis force.
"""
function assemble!(P::HDProblem{T,V}; threads=false, kwargs...) where {T,V}

    vbasis = P.vbasis.el
    
    nu = length(vbasis)
    nmat = nu
    TJ = promote_type(coefficienttype(vbasis[1][1]),eltype(P.cmat))

    vbasis = vptype{TJ}.(vbasis)
     
    cmat = convert.(TJ, P.cmat)

    # pfun! = threads ? projectforcet! : projectforce!
    if threads
        RHST = zeros(eltype(P.RHS),size(P.RHS)...)
        if !(isorthonormal(P.vbasis) && isorthonormal(P.bbasis))
            LHST = zeros(eltype(P.LHS),size(P.LHS)...)
            projectforcet!(view(LHST, 1:nu, 1:nu), cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t
            P.LHS = sparse(LHST)
            println("assemble LHS done!")
        else
            P.LHS = one(P.LHS)
        end

        projectforcet!(view(RHST, 1:nu, 1:nu), cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...) #Ω×u
        println("assemble 2/Le ∫ uᵢ⋅Ω×uⱼ dV done!")

        P.RHS = sparse(RHST)
    else
        if !(isorthonormal(P.vbasis) && isorthonormal(P.bbasis)) 
            projectforce!(view(P.LHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, inertial) #∂u/∂t
        else
            P.LHS = one(P.LHS)
        end

        projectforce!(view(P.RHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, coriolis, P.Ω) #Ω×u
    end
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
    TJ = promote_type(coefficienttype(vbasis[1][1]),coefficienttype(bbasis[1][1]),eltype(P.cmat))
    bbasis = vptype{TJ}.(bbasis)
    vbasis = vptype{TJ}.(vbasis)
     
    cmat = convert.(TJ, P.cmat)

    # pfun! = threads ? projectforcet! : projectforce!
    if threads
        RHST = zeros(eltype(P.RHS),size(P.RHS)...)
        if !(isorthonormal(P.vbasis) && isorthonormal(P.bbasis))
            LHST = zeros(eltype(P.LHS),size(P.LHS)...)

            projectforcet!(view(LHST, 1:nu, 1:nu), cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t
            if typeof(P.bbasis) <: InsulatingMFBasis
                ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
                LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
                ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(ls)))
                projectforcet_symmetric_neighbours!(view(LHST,nu+1:nmat,nu+1:nmat),cmat,bbasis,bbasis, inertial,LS,MS,ispt; kwargs...) #∂j/∂t
            else
                projectforce!(view(LHST, nu+1:nmat, nu+1:nmat), cmat, bbasis, bbasis, inertial; kwargs...) #∂j/∂t
            end
            P.LHS = sparse(LHST)
            println("assemble LHS done!")
        else
            P.LHS = one(P.LHS)
        end

        projectforcet!(view(RHST, 1:nu, 1:nu), cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...) #Ω×u
        println("assemble 2/Le ∫ uᵢ⋅Ω×uⱼ dV done!")
        projectforcet!(view(RHST, 1:nu, nu+1:nmat), cmat, vbasis, bbasis, lorentz, P.B0; kwargs...) #j×b
        println("assemble ∫ uᵢ⋅(∇×B₀×Bⱼ + ∇×Bⱼ×B₀) dV done!")
        projectforcet!(view(RHST, nu+1:nmat, 1:nu), cmat, bbasis, vbasis, advection, P.B0; kwargs...)
        println("assemble ∫ uᵢ⋅∇×uⱼ×B₀ done!")

        if !isinf(P.Lu)
            projectforcet!(view(RHST, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, b->1/P.Lu*diffusion(b); kwargs...)
            println("assemble 1/Lu ∫ Bᵢ⋅ΔBⱼ² dV done!")
        end

        P.RHS = sparse(RHST)
    else
        if !(isorthonormal(P.vbasis) && isorthonormal(P.bbasis)) 
            projectforce!(view(P.LHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, inertial) #∂u/∂t
            projectforce!(view(P.LHS, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, inertial) #∂j/∂t
        else
            P.LHS = one(P.LHS)
        end

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
