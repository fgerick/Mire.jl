### Functions to create matrices and assemble the full system.

abstract type MireProblem{T,V} end


"""
    HDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T, Vol}

Defines hydrodynamic problem.

Example:

```
N = 5
Ω = [0,0,1.0]
V = Ellipsoid(1.1,1.0,0.9)
problem = HDProblem(N,V,Ω,LebovitzBasis)
```
"""
mutable struct HDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T, Vol}
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
mutable struct MHDProblem{T,Vol<:Volume{T}} <: MireProblem{T, Vol}
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
    
    (T == Float64) && (N > 15) && @warn("N = $(N) with 64-bit floating point numbers will lead to inaccuricies! Consider using more accurate floats.")

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

    return MHDProblem(N, V, Ω, Le, T(Inf), B0, VB, BB; kwargs...)
end

"""
    assemble!(P::HDProblem{T,V}; threads=false, verbose=false, kwargs...) where {T,V}

Assembles the matrices `P.LHS` and `P.RHS`, i.e. projecting the velocity basis
`P.vbasis` on the inertial acceleration and Coriolis force.
"""
function assemble!(P::HDProblem{T,V}; threads=false, verbose=false, kwargs...) where {T,V}

    if threads
        assemblet!(P; verbose, kwargs...)
    else
        assembles!(P; verbose, kwargs...)
    end

    return nothing
end

function assemblet!(P::HDProblem{T,V}; verbose=false, kwargs...) where {T,V}

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

    return nothing
end

function assembles!(P::HDProblem{T,V}; verbose=false, kwargs...) where {T,V} 

    vbasis = P.vbasis.el

    if !isorthonormal(P.vbasis) 
        P.LHS = projectforce(vbasis, vbasis, P.cmat, inertial) #∂u/∂t
    else
        P.LHS = one(P.LHS)
    end

    P.RHS = projectforce(vbasis, vbasis, P.cmat, coriolis, P.Ω) #Ω×u

    return nothing
end

"""
    assemble!(P::MHDProblem{T,V}) where {T,V}

Assembles the matrices `P.LHS` and `P.RHS`, i.e. projecting the velocity and
magnetic field bases on the inertial acceleration, Coriolis force, Lorentz force
and mgnetic advection.
"""
function assemble!(P::MHDProblem{T,V}; threads=false, verbose=false, kwargs...) where {T,V}
    if threads
        assemblet!(P; verbose, kwargs...)
    else
        assembles!(P; verbose, kwargs...)
    end
    return nothing
end

function assemblet!(P::MHDProblem{T,V}; verbose=false, kwargs...) where {T,V}
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu + nb
    cmat = P.cmat

    nt = Threads.nthreads()
    
    if !(Mire.isorthonormal(P.vbasis) && Mire.isorthonormal(P.bbasis))

        itemps = [Int[] for i=1:nt]
        jtemps = [Int[] for i=1:nt]
        valtemps = [eltype(P.LHS)[] for i=1:nt]
        @sync begin 
            projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t
            
            if typeof(P.bbasis) <: Union{InsulatingMFBasis, InsMFONBasis, InsMFONCBasis, InsulatingMFCBasis}
                ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
                LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
                ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(lstor)))
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
            if typeof(P.bbasis) <: Union{InsulatingMFBasis, InsMFONBasis, InsMFONCBasis, InsulatingMFCBasis}
                ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
                LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
                ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(lstor)))
                Mire.projectforcet_symmetric_neighbours!(nu,nu,itemps,jtemps,valtemps,cmat,bbasis,bbasis, b->1/P.Lu*diffusion(b),LS,MS,ispt) #∂j/∂t
            else
                projectforcet!(nu, nu, itemps, jtemps, valtemps,  cmat, bbasis, bbasis, b->1/P.Lu*diffusion(b); kwargs...) #∂j/∂t
            end
            verbose && println("assemble 1/Lu ∫ Bᵢ⋅ΔBⱼ² dV done!")
        end
    @sync true
    P.RHS = sparse(vcat(itemps...),vcat(jtemps...),vcat(valtemps...), nmat, nmat)

end

function assembles!(P::MHDProblem{T,V}; kwargs...) where {T,V}
    
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    nu = length(vbasis)
    nmat = size(P.LHS,1)

    itemps,jtemps,valtemps = Int[], Int[], eltype(P.LHS)[]
    if !(isorthonormal(P.vbasis) && isorthonormal(P.bbasis)) 
        projectforce!(0, 0, itemps, jtemps, valtemps, P.cmat, vbasis, vbasis, inertial; kwargs...) #∂u/∂t

        if typeof(P.bbasis) <: Union{InsulatingMFBasis, InsMFONBasis, InsMFONCBasis, InsulatingMFCBasis}
            ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
            LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
            ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(lstor)))
            Mire.projectforce_symmetric_neighbours!(nu,nu,itemps,jtemps,valtemps,P.cmat,bbasis,bbasis, inertial,LS,MS,ispt) #∂j/∂t
        else
            projectforce!(nu, nu, itemps, jtemps, valtemps, P.cmat, bbasis, bbasis, inertial; kwargs...) #∂j/∂t
        end
        P.LHS = sparse(itemps,jtemps,valtemps, nmat, nmat)
    else
        P.LHS = one(P.LHS)
    end

    itemps, jtemps, valtemps = Int[], Int[], eltype(P.RHS)[]
    projectforce!(0,0, itemps, jtemps, valtemps, P.cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...) #Ω×u
    projectforce!(0, nu, itemps, jtemps, valtemps, P.cmat, vbasis, bbasis, lorentz, P.B0; kwargs...) #j×b
    projectforce!(nu, 0, itemps, jtemps, valtemps, P.cmat,  bbasis, vbasis, advection, P.B0; kwargs...)     

    if !isinf(P.Lu)
        projectforce!(nu, nu, itemps, jtemps, valtemps, P.cmat, bbasis, bbasis, b->1/P.Lu*diffusion(b); kwargs...)
    end

    P.RHS = sparse(itemps, jtemps, valtemps, nmat, nmat)

    return nothing
end

function assemblespecialized!(P::MHDProblem{T,V}, lb0, mb0, b0isp; verbose=false, kwargs...) where {T,V}
    @assert typeof(P.vbasis) <: Union{QGIMBasis,QGRIMBasis}
    @assert typeof(P.bbasis) <: Union{InsMFONCBasis, InsulatingMFBasis, InsulatingMFCBasis}
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu + nb
    cmat = P.cmat
    

    P.LHS = one(P.LHS) #only true for the given (normalized) bases
    
    #right hand side:
    ls,ms,ns,lstor,mstor,nstor = Mire.LMN(P.bbasis)

    N = P.N
    msu = Mire.NM(P.vbasis)[2]
    

    LS,MS = vcat(ls,lstor), vcat(ms, mstor)
    ISPT = vcat(fill(true,length(ls)),fill(false,length(lstor)))

    nt = Threads.nthreads()
    itemps = [Int[] for i=1:nt]
    jtemps = [Int[] for i=1:nt]
    valtemps = [eltype(P.LHS)[] for i=1:nt]
    @sync begin
    if typeof(P.vbasis) <: QGIMBasis
        projectcoriolisqgt!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, P.Ω; kwargs...) #Ω×u
    else
        projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vbasis, vbasis, coriolis, P.Ω; kwargs...)
    end
    verbose && println("assemble 2/Le ∫ uᵢ⋅Ω×uⱼ dV done!")

    projectlorentzqgt!(0, nu, itemps, jtemps, valtemps, cmat, vbasis, bbasis, P.B0, LS, MS, msu, lb0, mb0, ISPT, b0isp; kwargs...) #j×b
    verbose && println("assemble ∫ uᵢ⋅(∇×B₀×Bⱼ + ∇×Bⱼ×B₀) dV done!")

    projectinductionqgt!(nu, 0, itemps, jtemps, valtemps, cmat, bbasis, vbasis, P.B0, LS, MS, msu, lb0, mb0, ISPT, b0isp; kwargs...)
    verbose && println("assemble ∫ Bᵢ⋅∇×uⱼ×B₀ done!")

    if !isinf(P.Lu)
        ls,ms,ns,lstor,mstor,nstor = LMN(P.bbasis)
        LS,MS,NS = vcat(ls,lstor), vcat(ms,mstor), vcat(ns,nstor)
        ispt = vcat(zeros(Bool,length(ls)),ones(Bool,length(lstor)))
        Mire.projectforcet_symmetric_neighbours!(nu,nu,itemps,jtemps,valtemps,cmat,bbasis,bbasis, b->1/P.Lu*diffusion(b),LS,MS,ispt) #∂j/∂t
        verbose && println("assemble 1/Lu ∫ Bᵢ⋅ΔBⱼ² dV done!")
    end
    end
    P.RHS = sparse(vcat(itemps...),vcat(jtemps...),vcat(valtemps...), nmat, nmat)

    return nothing
end



