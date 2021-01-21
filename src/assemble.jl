### Functions to create matrices and assemble the full system.

abstract type MireProblem{T,V} end



mutable struct HDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T,Vol}
    N::Int
    V::Vol
    Ω::Vector{T}
    vbasis::VectorBasis{T,Vol}
    cmat::Array{T,3}
    LHS::AbstractMatrix{T}
    RHS::AbstractMatrix{T}
end

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

function HDProblem(N::Int, V::Volume{T}, Ω::Vector{T}, ::Type{VB}) where {T <: Number, VB <: VectorBasis}
    vbasis = VB(N,V)
    cmat = cacheint(N,V)
    n = length(vbasis.el)
    LHS = spzeros(T,n,n)
    RHS = spzeros(T,n,n)
     
    return HDProblem(N, V, Ω, vbasis, cmat, LHS, RHS)
end

function MHDProblem(N::Int, V::Volume{T}, Ω::Vector{T}, B0, ::Type{VB}, ::Type{BB}) where {T <: Number, VB <: VectorBasis, BB <: VectorBasis}
    vbasis = VB(N,V)
    bbasis = BB(N,V)
    cmat = cacheint(N,V)
    nu = length(vbasis.el)
    nb = length(bbasis.el)
    n = nu + nb
    LHS = spzeros(T,n,n)
    RHS = spzeros(T,n,n)
     
    return MHDProblem(N, V, Ω, vptype{T}(B0), vbasis, bbasis, cmat, LHS, RHS)
end


function assemble!(P::HDProblem{T,V}) where {T,V}
    
    projectforce!(P.LHS,P.cmat,P.vbasis.el,inertial)
    projectforce!(P.RHS,P.cmat,P.vbasis.el,coriolis,P.Ω)

    return nothing
end


function assemble!(P::MHDProblem{T,V}) where {T,V}
    vbasis = P.vbasis.el
    bbasis = P.bbasis.el
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu+nb
    projectforce!(view(P.LHS,1:nu,1:nu),P.cmat,vbasis,vbasis, inertial) #∂u/∂t
    projectforce!(view(P.LHS,nu+1:nmat,nu+1:nmat),P.cmat,bbasis,bbasis, inertial) #∂j/∂t
    projectforce!(view(P.RHS,1:nu,1:nu),P.cmat,vbasis,vbasis,coriolis,P.Ω) #Ω×u
    projectforce!(view(P.RHS,1:nu,nu+1:nmat),P.cmat,vbasis,bbasis,lorentz,P.B0) #j×b
    projectforce!(view(P.RHS,nu+1:nmat,1:nu),P.cmat,bbasis,vbasis,advection,P.B0)
    nothing
end

"""
    assemblemhd(N::Int, a::T, b::T, c::T, Ω, b0; dtype::DataType=BigFloat) where T

Assemble the sparse matrices of the MHD mode problem. Returns right hand side `A`,
left hand side `B` and basis vectors `vbasis`.

#Arguments:
- `N`: maximum monomial degree
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `b0`: mean magnetic field vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
"""
function assemblemhd(N::Int,a::T,b::T,c::T,Ω,b0;
                     cmat = cacheint(N,a,b,c)) where T

    vbasis = vel(N,a,b,c)
    bbasis = vbasis
    A,B = assemblemhd(Ω,b0,vbasis,bbasis,cmat)
    return A,B, vbasis
end

function assemblemhd(Ω,b0,vbasis,bbasis,cmat::Array{T,3}) where T
    nu = length(vbasis)
    nb = length(bbasis)
    nmat = nu+nb
    A = spzeros(T,nmat,nmat)
    B = spzeros(T,nmat,nmat)
    assemblemhd!(A,B,cmat,vbasis,bbasis,Ω,b0)
    return A,B
end



"""
    assemblehd_hybrid(N2D::Int, N3D::Int, a::T, b::T, c::T, Ω ; dtype::DataType=BigFloat) where T

Assemble the sparse matrices of the hybrid QG and 3D MHD mode problem.
Returns right hand side `A`,left hand side `B` and basis vectors `bbasis` and `vbasis_qg`.

#Arguments:
- `N2D`: maximum monomial degree of QG velocity
- `N3D`: maximum monomial degree of 3D magnetic field
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
"""
function assemblemhd_hybrid(N2D::Int,N3D::Int,a::T,b::T,c::T,Ω,b0;
                     cmat = cacheint(N3D,a,b,c)) where T
    # n_mat = Mire.n_u(N3D)
    bbasis = Mire.vel(N3D,a,b,c)
    n_mat = length(bbasis)
    vbasis_qg = Mire.qg_vel(N2D,a,b,c)
    n_mat_qg = length(vs_qg)


    nmat=n_mat+n_mat_qg
    A = spzeros(T,nmat,nmat)
    B = spzeros(T,nmat,nmat)

    assemblemhd!(A,B,cmat,vbasis_qg,bbasis,Ω,b0)
    return A,B, bbasis, vbasis_qg
end

"""
    assemblehd_qg(N2D::Int, a::T, b::T, c::T, Ω ; dtype::DataType=BigFloat) where T

Assemble the sparse matrices of the QG MHD mode problem.
Returns right hand side `A`,left hand side `B` and basis vectors `vs_qg`.

#Arguments:
- `N2D`: maximum monomial degree of QG velocity/magnetic field
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
"""
function assemblemhd_qg(N2D::Int, a::T, b::T, c::T, Ω, b0;
                     cmat = cacheint(N2D,a,b,c)) where T
    vs_qg = Mire.qg_vel(N2D, a, b, c)
    n_mat_qg = length(vs_qg)

    nmat = 2n_mat_qg
    A = spzeros(T, nmat, nmat)
    B = spzeros(T, nmat, nmat)
    assemblemhd!(A,B,cmat,vs_qg,vs_qg,Ω,b0)

    return A, B, vs_qg
end


#Quagmire (2D reduced equations)

function assemblemhd_quag(N::Int, a::T, b::T, c::T, Ω::T, A0;
                     cmat = cacheint2D(N,a,b,c)) where T
    vs_qg = Mire.qg_vel(N, a, b, c)
    n_mat_qg = length(vs_qg)

    nmat = 2n_mat_qg
    A = spzeros(T, nmat, nmat)
    B = spzeros(T, nmat, nmat)
    projectforce_2D!(view(A, 1:n_mat_qg, 1:n_mat_qg),           N, cmat, Mire.poly_inertial,a,b,c)
    projectforce_2D!(view(A, n_mat_qg+1:nmat, n_mat_qg+1:nmat), N, cmat, Mire.poly_inertialmag,a,b,c)
    projectforce_2D!(view(B, 1:n_mat_qg, 1:n_mat_qg),           N, cmat, Mire.poly_coriolis,a,b,c, Ω)
    projectforce_2D!(view(B, 1:n_mat_qg, n_mat_qg+1:nmat),      N, cmat, Mire.poly_lorentz,a,b,c, A0)
    projectforce_2D!(view(B, n_mat_qg+1:nmat, 1:n_mat_qg),      N, cmat, Mire.poly_advection,a,b,c, A0)

    return A, B, vs_qg
end



function projectforce_2D!(A::AbstractArray{T,2},N::Integer, cmat, polyfun::Function,a::T,b::T,c::T, args...) where T

    combos = qg_combos(N)

    n_A = length(combos)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A

    @inbounds for j=1:n_A
        p = polyfun(combos[j]...,a,b,c,args...)
        for i=1:n_A
            A[i,j] = inner_product_2D(cmat,p,Π(combos[i]...),a,b,c)
        end
    end
end


function projectforce_2D(N::Integer,cmat, polyfun::Function, a::T,b::T,c::T, args...) where T
    n_combos = n_unknown(N)
    A = spzeros(T,n_combos,n_combos)
    projectforce_2D!(A,N,cmat, polyfun,a,b,c,args...)
    return A
end
