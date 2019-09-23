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

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A
    @assert length(vs)==n_A


    @inbounds for j=1:n_A
        f = forcefun(vs[j],args...) #calculate f(uⱼ)
        for i=1:n_A
            A[i,j] = inner_product(cmat,vs[i],f; kwargs...)
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
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    projectforce!(A,cmat,vs,N ,forcefun,args...;kwargs...)
    return A
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
                     dtype::DataType=BigFloat, kwargs...) where T
    n_mat = n_u(N)
    vs = vel(N,a,b,c)
    cmat = cacheint(N,a,b,c; dtype=dtype)

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
function assemblehd(N::Int,a::T,b::T,c::T,Ω ;
                    dtype::DataType=BigFloat, kwargs...) where T
    cmat = cacheint(N,a,b,c; dtype=dtype)
    vs = vel(N,a,b,c)

    A = projectforce(N,cmat,vs,inertial; kwargs...)
    B = projectforce(N,cmat,vs,coriolis,Ω; kwargs...)
    return A,B, vs
end
