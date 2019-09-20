### Functions to create matrices and assemble the full system.

"""
    projectforce!(A::AbstractArray{T, 2}, cmat::Array{T, 3}, vs::Array{Array{P, 1}, 1}, N::Integer, forcefun::Function, a::T, b::T, c::T, args...; kwargs...) where {T, P <: Polynomial{T}}

DOCSTRING

#Arguments:
- `A`: DESCRIPTION
- `cmat`: DESCRIPTION
- `vs`: DESCRIPTION
- `N`: DESCRIPTION
- `forcefun`: DESCRIPTION
- `a`: DESCRIPTION
- `b`: DESCRIPTION
- `c`: DESCRIPTION
- `args`: DESCRIPTION
- `kwargs`: DESCRIPTION
"""
function projectforce!(A::AbstractArray{T,2},cmat::Array{T,3},vs::Array{Array{P,1},1},
            N::Integer, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T, P <: Polynomial{T}}

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A
    @assert length(vs)==n_A


    @inbounds for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            A[i,j] = inner_product(cmat,vs[i],f; kwargs...)
        end
    end
end


"""
    projectforce(N::Integer, vs, forcefun::Function, a::T, b::T, c::T, args...) where T

Allocates new matrix `A` and fills elements by calling
projectforce!(A,vs,N,forcefun,a,b,c, args...).

Cached version:
projectforce(N,cmat,vs,forcefun,a,b,c, args...)

where `cmat[i,j,k]` contains the integrals of monomials xⁱyʲzᵏ.

#Arguments:
- `N`: DESCRIPTION
- `vs`: DESCRIPTION
- `forcefun`: DESCRIPTION
- `a`: DESCRIPTION
- `b`: DESCRIPTION
- `c`: DESCRIPTION
- `args`: DESCRIPTION
"""
function projectforce(N::Integer,vs, forcefun::Function,a::T,b::T,c::T, args...) where T
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    projectforce!(A,vs,N ,forcefun,a,b,c,args...)
    return A
end

function projectforce(N::Integer,cmat::Array{T,3},vs::Array{Array{P,1},1}, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T, P <: Polynomial{T}}
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    projectforce!(A,cmat,vs,N ,forcefun,a,b,c,args...;kwargs...)
    return A
end

"""
    assemblemhd(N::Int, cmat::Array{T, 3}, a::T, b::T, c::T, Ω, b0; kwargs...) where T

DOCSTRING

#Arguments:
- `N`: DESCRIPTION
- `cmat`: DESCRIPTION
- `a`: DESCRIPTION
- `b`: DESCRIPTION
- `c`: DESCRIPTION
- `Ω`: DESCRIPTION
- `b0`: DESCRIPTION
- `kwargs`: DESCRIPTION
"""
function assemblemhd(N::Int,cmat::Array{T,3},a::T,b::T,c::T,Ω,b0; kwargs...) where T
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    A[1:n_mat,1:n_mat] .= projectforce(N,cmat,vs,inertial,a,b,c; kwargs...)
    A[n_mat+1:end,n_mat+1:end] .= projectforce(N,cmat,vs,inertialmag,a,b,c; kwargs...)

    B[1:n_mat,1:n_mat] .= projectforce(N,cmat,vs,coriolis,a,b,c,Ω; kwargs...)
    B[1:n_mat,n_mat+1:end] .= projectforce(N,cmat,vs,lorentz,a,b,c,b0; kwargs...)

    B[n_mat+1:end,1:n_mat] .= projectforce(N,cmat,vs,advection,a,b,c,b0; kwargs...)

    return A,B, vs
end

using AutomaticDocstrings

"""
    assemblehd(N::Int, cmat, a::T, b::T, c::T, Ω) where T

DOCSTRING

#Arguments:
- `N`: DESCRIPTION
- `cmat`: DESCRIPTION
- `a`: DESCRIPTION
- `b`: DESCRIPTION
- `c`: DESCRIPTION
- `Ω`: DESCRIPTION
"""
function assemblehd(N::Int,cmat,a::T,b::T,c::T,Ω) where T

    vs = vel(N,a,b,c)

    A = projectforce(N,cmat,vs,inertial,a,b,c)
    B = projectforce(N,cmat,vs,coriolis,a,b,c,Ω)
    return A,B, vs
end
