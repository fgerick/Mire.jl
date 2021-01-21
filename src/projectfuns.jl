## low level functions that project the basis elements on a force function (forcefun)

"""
projectforce!(A::AbstractArray{T, 2}, cmat::Array{T, 3}, vs::Array{Array{P, 1}, 1}, N::Integer, forcefun::Function, a::T, b::T, c::T, args...) where {T, P <: Polynomial{T}}

DOCSTRING

#Arguments:
- `A`: pre-allocated array
- `cmat`: pre-cached monomial integration values
- `vs`: basis vectors
- `N`: maximum monomial degree
- `forcefun`: function of the force, e.g. coriolis
- `args`: other arguments needed for `forcefun`
"""
function projectforce!(A::AbstractArray{T,2},cmat::Array{T,3},vs::Array{Array{P,1},1}, 
                        forcefun::Function, args...) where {T, P <: Polynomial{T}}
    projectforce!(A, cmat, vs, vs, forcefun, args...)
end

function projectforce!(A::AbstractArray{T,2},cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},
        forcefun::Function, args...) where {T, P <: Polynomial{T}}

n_1 = length(vs_i)
n_2 = length(vs_j)
@assert n_1 == size(A,1)
@assert n_2 == size(A,2)

@inbounds for j=1:n_2
    f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
    for i=1:n_1
        A[i,j] = inner_product_real(cmat,vs_i[i],f)
    end
end
end


"""
projectforce(N::Integer,cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},
forcefun::Function,a::T,b::T,c::T, args...) where {T, P <: Polynomial{T}}

Allocates new matrix `A` and fills elements by calling
projectforce!(A,cmat,vs_i,vs_j,forcefun, args...)

where `cmat[i,j,k]` contains the integrals of monomials xⁱyʲzᵏ.

#Arguments:
- `N`: maximum monomial degree
- `vs_i`: basis vectors to project on
- `vs_j`: basis vectors used for `forcefun`
- `forcefun`: function of the force, e.g. coriolis
- `args`: other arguments needed for `forcefun`
"""
function projectforce(cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},
    forcefun::Function, args...) where {T, P <: Polynomial{T}}

n_1 = length(vs_i)
n_2 = length(vs_j)

A = spzeros(T,n_1,n_2)

projectforce!(A, cmat, vs_i, vs_j, forcefun, args...)
return A
end

projectforce(cmat::Array{T,3},vs::Array{Array{P,1},1},forcefun::Function, args...) where {T, P <: Polynomial{T}} = projectforce(cmat,vs,vs,forcefun,args...)

