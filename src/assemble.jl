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
mutable struct MHDProblem{T,Vol<:Volume{T}} <: MireProblem{T,Vol}
    N::Int
    V::Vol
    Ω::Vector{T}
    Le::T
    Lu::T
    B0 #::Union{vptype{T},vptype{Complex{T}}}
    vbasis #::VectorBasis{T,Vol}
    bbasis #::VectorBasis{T,Vol}
    cmat #::Array{T,3}
    LHS #::Union{AbstractMatrix{T},AbstractMatrix{Complex{T}}}
    RHS #::Union{AbstractMatrix{T},AbstractMatrix{Complex{T}}}
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
                    ::Type{BB};
                    kwargs...
                ) where {T<:Number,VB<:VectorBasis,BB<:VectorBasis}

    vbasis = VB(N, V)
    bbasis = BB(N, V; kwargs...)
    cmat = cacheint(N, V)
    nu = length(vbasis.el)
    nb = length(bbasis.el)
    n = nu + nb
    TM = promote_type(coefficienttype(vbasis.el[1][1]),coefficienttype(bbasis.el[1][1]))
    LHS = spzeros(TM, n, n)
    RHS = spzeros(TM, n, n)
    if norm(Ω) != 1
        Ω /= Le*norm(Ω)
    else
        Ω = Ω*1/Le
    end
    B01=vptype{TM}(B0)
    return MHDProblem(N, V, Ω, Le, Lu, B01, vbasis, bbasis, cmat, LHS, RHS)
end

function MHDProblem(
    N::Int,
    V::Volume{T},
    Ω::Vector{T},
    Le::T,
    B0,
    ::Type{VB},
    ::Type{BB};
    kwargs...
) where {T<:Number,VB<:VectorBasis,BB<:VectorBasis}

    # vbasis = VB(N, V)
    # bbasis = BB(N, V; kwargs...)
    # cmat = cacheint(N, V)
    # nu = length(vbasis.el)
    # nb = length(bbasis.el)
    # n = nu + nb
    # TM = promote_type(coefficienttype(vbasis.el[1][1]),coefficienttype(bbasis.el[1][1]))
    # @show TM
    # LHS = spzeros(TM, n, n)
    # RHS = spzeros(TM, n, n)

    return MHDProblem(N, V, Ω, Le, T(Inf), B0, VB, BB; kwargs...)
end

"""
    assemble!(P::HDProblem{T,V}) where {T,V}

Assembles the matrices `P.LHS` and `P.RHS`, i.e. projecting the velocity basis
`P.vbasis` on the inertial acceleration and Coriolis force.
"""
function assemble!(P::HDProblem{T,V}; threads=false, verbose=false, kwargs...) where {T,V}

    vbasis = P.vbasis.el
    
    nu = length(vbasis)
    nmat = nu
    TJ = promote_type(coefficienttype(vbasis[1][1]),eltype(P.cmat))

    vbasis = vptype{TJ}.(vbasis)
     
    cmat = convert.(TJ, P.cmat)

    # pfun! = threads ? projectforcet! : projectforce!
    if threads
        nt = Threads.nthreads()
        #LHS
        if !isorthonormal(P.vbasis)
            itemps = [Int[] for i=1:nt]
            jtemps = [Int[] for i=1:nt]
            valtemps = [eltype(P.LHS)[] for i=1:nt]
    
            @sync projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t 
            P.LHS = sparse(vcat(itemps...),vcat(jtemps...),vcat(valtemps...), nu, nu)
            verbose && println("assemble LHS done!")
        else
            P.LHS = one(P.LHS)
        end
        
        
        #RHS
        itemps = [Int[] for i=1:nt]
        jtemps = [Int[] for i=1:nt]
        valtemps = [eltype(P.LHS)[] for i=1:nt]

        @sync projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...) #Ω×u
        P.RHS = sparse(vcat(itemps...),vcat(jtemps...),vcat(valtemps...), nu, nu)
        verbose && println("assemble 2/Le ∫ uᵢ⋅Ω×uⱼ dV done!")

    else
        if !isorthonormal(P.vbasis) 
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
function assemble!(P::MHDProblem{T,V}; threads=false, verbose=false, kwargs...) where {T,V}
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu + nb
    cmat = P.cmat

    
    if threads
        nt = Threads.nthreads()

        
        if !(Mire.isorthonormal(P.vbasis) && Mire.isorthonormal(P.bbasis))
            # LHST = view(P.LHS, 1:nu, 1:nu) #zeros(eltype(P.LHS),nu,nu)

            itemps = [Int[] for i=1:nt]
            jtemps = [Int[] for i=1:nt]
            valtemps = [eltype(P.LHS)[] for i=1:nt]
            @sync begin 
                projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t
                
                if typeof(P.bbasis) <: Union{InsulatingMFBasis, InsMFONBasis, InsMFONCBasis, InsMFCBasis}
                    ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
                    LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
                    ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(ls)))
                    projectforcet_symmetric_neighbours!(nu,nu,itemps,jtemps,valtemps,cmat,bbasis,bbasis, inertial,LS,MS,ispt) #∂j/∂t
                else
                    projectforcet!(nu, nu, itemps, jtemps, valtemps, cmat, bbasis, bbasis, inertial; kwargs...) #∂j/∂t
                end
                verbose && println("assemble LHS done!")
            end
                P.LHS = sparse(vcat(itemps...),vcat(jtemps...),vcat(valtemps...), nmat, nmat)
        else
            P.LHS = one(P.LHS)
        end

        #right hand side:

        itemps = [Int[] for i=1:nt]
        jtemps = [Int[] for i=1:nt]
        valtemps = [eltype(P.LHS)[] for i=1:nt]

            projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...) #Ω×u
            verbose && println("assemble 2/Le ∫ uᵢ⋅Ω×uⱼ dV done!")

            projectforcet!(0, nu, itemps, jtemps, valtemps, cmat, vbasis, bbasis, lorentz, P.B0; kwargs...) #j×b
            verbose && println("assemble ∫ uᵢ⋅(∇×B₀×Bⱼ + ∇×Bⱼ×B₀) dV done!")

            projectforcet!(nu, 0, itemps, jtemps, valtemps, cmat, bbasis, vbasis, advection, P.B0; kwargs...)
            verbose && println("assemble ∫ Bᵢ⋅∇×uⱼ×B₀ done!")

            if !isinf(P.Lu)
                if typeof(P.bbasis) <: Union{InsulatingMFBasis, InsMFONBasis, InsMFONCBasis, InsMFCBasis}
                    ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
                    LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
                    ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(ls)))
                    Mire.projectforcet_symmetric_neighbourst!(nu,nu,itemps,jtemps,valtemps,cmat,bbasis,bbasis, b->1/P.Lu*diffusion(b),LS,MS,ispt) #∂j/∂t
                else
                    projectforcet!(nu, nu, itemps, jtemps, valtemps,  cmat, bbasis, bbasis, b->1/P.Lu*diffusion(b); kwargs...) #∂j/∂t
                end
                verbose && println("assemble 1/Lu ∫ Bᵢ⋅ΔBⱼ² dV done!")
            end
        @sync true
        P.RHS = sparse(vcat(itemps...),vcat(jtemps...),vcat(valtemps...), nmat, nmat)
    else #no threads version
        if !(isorthonormal(P.vbasis) && isorthonormal(P.bbasis)) 
            projectforce!(view(P.LHS, 1:nu, 1:nu), P.cmat, vbasis, vbasis, inertial) #∂u/∂t

            if typeof(P.bbasis) <: Union{InsulatingMFBasis, InsMFONBasis, InsMFONCBasis, InsMFCBasis}
                ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
                LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
                ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(ls)))
                Mire.projectforce_symmetric_neighbours!(view(P.LHS,nu+1:nmat,nu+1:nmat),cmat,bbasis,bbasis, inertial,LS,MS,ispt; kwargs...) #∂j/∂t
            else
                projectforce!(view(P.LHS, nu+1:nmat, nu+1:nmat), P.cmat, bbasis, bbasis, inertial; kwargs...) #∂j/∂t
            end
            
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

    return nothing
end


function assemblerhs!(P::MHDProblem{T,V}; threads=false, kwargs...) where {T,V}

end

function normalizebasis!(P::MireProblem{T}; n_cache::Int = 10^6) where T
    ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache)
    el = map(u->u/sqrt(Mire.inner_product!(ptemp,u,u,P.cmat)) , P.vbasis.el)
    P.vbasis = typeof(P.vbasis)(P.N,P.V,el)

    if typeof(P) <: MHDProblem

       normfac = map(b->sqrt(Mire.inner_product!(ptemp,b,b,P.cmat)) , P.bbasis.el)
       el = map(b->b/normfac,P.bbasis.el)
       P.bbasis = typeof(P.bbasis)(P.N,P.V,el)
    end
    return normfac
end
