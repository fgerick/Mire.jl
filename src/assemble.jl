### Functions to create matrices and assemble the full system.

"""
    projectforce!(A::AbstractArray{T, 2}, cmat::Array{T, 3}, vs::Array{Array{P, 1}, 1}, N::Integer, forcefun::Function, a::T, b::T, c::T, args...; kwargs...) where {T, P <: Polynomial{T}}

DOCSTRING

#Arguments:
- `A`: pre-allocated array
- `cmat`: pre-cached monomial integration values
- `vs`: basis vectors
- `N`: maximum monomial degree
- `forcefun`: function of the force, e.g. coriolis
- `args`: other arguments needed for `forcefun`
- `kwargs`: other keyword arguments
"""
function projectforce!(A::AbstractArray{T,2},cmat::Array{T,3},vs::Array{Array{P,1},1},
            N::Integer, forcefun::Function, args...; kwargs...) where {T, P <: Polynomial{T}}
    projectforce!(A, cmat, vs, vs, forcefun, args...; kwargs...)
end

function projectforce!(A::AbstractArray{T,2},cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},
            forcefun::Function, args...; kwargs...) where {T, P <: Polynomial{T}}

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    @assert n_1 == size(A,1)
    @assert n_2 == size(A,2)

    @inbounds for j=1:n_2
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        for i=1:n_1
            A[i,j] = Mire.inner_product(cmat,vs_i[i],f; kwargs...)
        end
    end
end

"""
    projectforce(N::Integer,cmat::Array{T,3},vs::Array{Array{P,1},1},
    forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T, P <: Polynomial{T}}

Allocates new matrix `A` and fills elements by calling
projectforce!(A,cmat,vs,forcefun, args...; kwargs...)

where `cmat[i,j,k]` contains the integrals of monomials xⁱyʲzᵏ.

#Arguments:
- `N`: maximum monomial degree
- `vs`: basis vectors
- `forcefun`: function of the force, e.g. coriolis
- `args`: other arguments needed for `forcefun`
"""
function projectforce(N::Integer,cmat::Array{T,3},vs::Array{Array{P,1},1},
        forcefun::Function, args...; kwargs...) where {T, P <: Polynomial{T}}
    n_combos = length(vs)
    A = spzeros(T,n_combos,n_combos)
    projectforce!(A, cmat, vs, N, forcefun, args...; kwargs...)
    return A
end

function projectforce(cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},
        forcefun::Function, args...; kwargs...) where {T, P <: Polynomial{T}}

    n_1 = length(vs_i)
    n_2 = length(vs_j)

    A = spzeros(T,n_1,n_2)

    projectforce!(A, cmat, vs_i, vs_j, forcefun, args...; kwargs...)
    return A
end


"""
    assemblehd(N::Int, a::T, b::T, c::T, Ω ; dtype::DataType=BigFloat, kwargs...) where T

Assemble the sparse matrices of the MHD mode problem. Returns right hand side `A`,
left hand side `B` and basis vectors `vs`.

#Arguments:
- `N`: maximum monomial degree
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
- `kwargs`: other keyword arguments passed to lower functions
"""
function assemblehd(N::Int,a::T,b::T,c::T,Ω ; kwargs...) where T
    cmat = cacheint(N,a,b,c)
    vs = vel(N,a,b,c)

    A = projectforce(N,cmat,vs,inertial; kwargs...)
    B = projectforce(N,cmat,vs,coriolis,Ω; kwargs...)
    return A,B, vs
end

"""
    assemblemhd(N::Int, a::T, b::T, c::T, Ω, b0; dtype::DataType=BigFloat, kwargs...) where T

Assemble the sparse matrices of the MHD mode problem. Returns right hand side `A`,
left hand side `B` and basis vectors `vs`.

#Arguments:
- `N`: maximum monomial degree
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `b0`: mean magnetic field vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
- `kwargs`: other keyword arguments passed to lower functions
"""
function assemblemhd(N::Int,a::T,b::T,c::T,Ω,b0;
                     cmat = cacheint(N,a,b,c), kwargs...) where T
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)
    projectforce!(view(A,1:n_mat,1:n_mat),cmat,vs,N, inertial; kwargs...)
    projectforce!(view(A,n_mat+1:2n_mat,n_mat+1:2n_mat),cmat,vs,N, inertialmag; kwargs...)
    projectforce!(view(B,1:n_mat,1:n_mat),cmat,vs,N,coriolis,Ω; kwargs...)
    projectforce!(view(B,1:n_mat,n_mat+1:2n_mat),cmat,vs,N,lorentz,b0; kwargs...)
    projectforce!(view(B,n_mat+1:2n_mat,1:n_mat),cmat,vs,N,advection,b0; kwargs...)

    return A,B, vs
end


"""
    assemblehd_hybrid(N2D::Int, N3D::Int, a::T, b::T, c::T, Ω ; dtype::DataType=BigFloat, kwargs...) where T

Assemble the sparse matrices of the hybrid QG and 3D MHD mode problem.
Returns right hand side `A`,left hand side `B` and basis vectors `vs` and `vs_qg`.

#Arguments:
- `N2D`: maximum monomial degree of QG velocity
- `N3D`: maximum monomial degree of 3D magnetic field
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
- `kwargs`: other keyword arguments passed to lower functions
"""
function assemblemhd_hybrid(N2D::Int,N3D::Int,a::T,b::T,c::T,Ω,b0;
                     cmat = cacheint(N3D,a,b,c), kwargs...) where T
    n_mat = Mire.n_u(N3D)
    vs = Mire.vel(N3D,a,b,c)
    vs_qg = Mire.qg_vel(N2D,a,b,c)
    n_mat_qg = length(vs_qg)


    nmat=n_mat+n_mat_qg
    A = spzeros(T,nmat,nmat)
    B = spzeros(T,nmat,nmat)
    projectforce!(view(A,1:n_mat_qg,1:n_mat_qg),            cmat, vs_qg, vs_qg, Mire.inertial; kwargs...)
    projectforce!(view(A,n_mat_qg+1:nmat,n_mat_qg+1:nmat),  cmat, vs, vs, Mire.inertialmag; kwargs...)
    projectforce!(view(B,1:n_mat_qg,1:n_mat_qg),            cmat, vs_qg, vs_qg, Mire.coriolis,Ω; kwargs...)
    projectforce!(view(B,1:n_mat_qg,n_mat_qg+1:nmat),       cmat, vs_qg, vs, Mire.lorentz,b0; kwargs...)
    projectforce!(view(B,n_mat_qg+1:nmat,1:n_mat_qg),       cmat, vs,vs_qg, Mire.advection,b0; kwargs...)

    return A,B, vs, vs_qg
end

"""
    assemblehd_qg(N2D::Int, a::T, b::T, c::T, Ω ; dtype::DataType=BigFloat, kwargs...) where T

Assemble the sparse matrices of the QG MHD mode problem.
Returns right hand side `A`,left hand side `B` and basis vectors `vs_qg`.

#Arguments:
- `N2D`: maximum monomial degree of QG velocity/magnetic field
- `a`: semi-axis x
- `b`: semi-axis y
- `c`: semi-axis z
- `Ω`: rotation vector
- `dtype`: datatype, default `BigFloat` for integration of monomials
- `kwargs`: other keyword arguments passed to lower functions
"""
function assemblemhd_qg(N2D::Int, a::T, b::T, c::T, Ω, b0;
                     cmat = cacheint(N2D,a,b,c), kwargs...) where T
    vs_qg = Mire.qg_vel(N2D, a, b, c)
    n_mat_qg = length(vs_qg)

    nmat = 2n_mat_qg
    A = spzeros(T, nmat, nmat)
    B = spzeros(T, nmat, nmat)
    projectforce!(view(A, 1:n_mat_qg, 1:n_mat_qg),           cmat, vs_qg, vs_qg, Mire.inertial; kwargs...)
    projectforce!(view(A, n_mat_qg+1:nmat, n_mat_qg+1:nmat), cmat, vs_qg, vs_qg, Mire.inertialmag; kwargs...)
    projectforce!(view(B, 1:n_mat_qg, 1:n_mat_qg),           cmat, vs_qg, vs_qg, Mire.coriolis, Ω; kwargs...)
    projectforce!(view(B, 1:n_mat_qg, n_mat_qg+1:nmat),      cmat, vs_qg, vs_qg, Mire.lorentz, b0; kwargs...)
    projectforce!(view(B, n_mat_qg+1:nmat, 1:n_mat_qg),      cmat, vs_qg, vs_qg, Mire.advection, b0; kwargs...)

    return A, B, vs_qg
end


#Quagmire

function assemblemhd_quag(N::Int, a::T, b::T, c::T, Ω::T, A0;
                     cmat = cacheint2D(N,a,b,c), kwargs...) where T
    vs_qg = Mire.qg_vel(N, a, b, c)
    n_mat_qg = length(vs_qg)

    nmat = 2n_mat_qg
    A = spzeros(T, nmat, nmat)
    B = spzeros(T, nmat, nmat)
    projectforce_2D!(view(A, 1:n_mat_qg, 1:n_mat_qg),           N, cmat, Mire.poly_inertial,a,b,c; kwargs...)
    projectforce_2D!(view(A, n_mat_qg+1:nmat, n_mat_qg+1:nmat), N, cmat, Mire.poly_inertialmag,a,b,c; kwargs...)
    projectforce_2D!(view(B, 1:n_mat_qg, 1:n_mat_qg),           N, cmat, Mire.poly_coriolis,a,b,c, Ω; kwargs...)
    projectforce_2D!(view(B, 1:n_mat_qg, n_mat_qg+1:nmat),      N, cmat, Mire.poly_lorentz,a,b,c, A0; kwargs...)
    projectforce_2D!(view(B, n_mat_qg+1:nmat, 1:n_mat_qg),      N, cmat, Mire.poly_advection,a,b,c, A0; kwargs...)

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
