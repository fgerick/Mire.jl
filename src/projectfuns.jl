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
function projectforce!(
    A::AbstractArray{T,2},
    cmat::Array{T,3},
    vs::Array{Array{P,1},1}, 
    forcefun::Function, 
    args...
    ) where {T, P <: Polynomial{T}}

    return projectforce!(A, cmat, vs, vs, forcefun, args...)
end

function projectforce!(
    A::AbstractArray{T,2},
    cmat::Array{T,3},
    vs_i::Array{Array{P,1},1},
    vs_j::Array{Array{P,1},1},
    forcefun::Function,
    args...
    ) where {T, P <: Polynomial{T}}

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

    return nothing
end


"""
    projectforce(N::Integer,cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},forcefun::Function,a::T,b::T,c::T, args...) where {T, P <: Polynomial{T}}

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
function projectforce(
    cmat::Array{T,3},
    vs_i::Array{Array{P,1},1},
    vs_j::Array{Array{P,1},1},
    forcefun::Function, 
    args...
    ) where {T, P <: Polynomial{T}}

    n_1 = length(vs_i)
    n_2 = length(vs_j)

    A = spzeros(T,n_1,n_2)

    projectforce!(A, cmat, vs_i, vs_j, forcefun, args...)

    return A
end

function projectforce(
    cmat::Array{T,3},
    vs::Array{Array{P,1},1},
    forcefun::Function, 
    args...
    ) where {T, P <: Polynomial{T}}

    return projectforce(cmat,vs,vs,forcefun,args...)
end



#vector type definition/abbreviation
ptype{T} = Polynomial{T,Term{T,Monomial{(x, y, z),3}},Array{Term{T,Monomial{(x, y, z),3}},1}}
vptype{T} = Vector{ptype{T}}

####
# inner product without simplifying the polynomial/merging polynomial terms.
# this way some overhead can be avoided, at the cost of increasingly large vectors
# of monomial terms.

function _mul_to_terms!(
    iter::Int,
    p::Vector{Term{T,Monomial{(x, y, z),3}}},
    u::ptype{T},
    v::ptype{T}
    ) where T

    t1 = terms(u)
    t2 = terms(v)
    @assert length(t1)*length(t2) < length(p)-iter "Pre-cached term array not large enough! $(length(t1)*length(t2)) >= $(length(p)-iter)"
    i = iter #number of used terms from preallocated terms-vector p
    @inbounds for ti in t1
        for tj in t2
            p[iter] = ti*tj
            iter += 1
        end
    end
    return iter
end

function dotpp!(
    p::Vector{Term{T,Monomial{(x, y, z),3}}},
    u::vptype{T},
    v::vptype{T}
    ) where T

    k=1
    @inbounds for i in 1:3
        k = _mul_to_terms!(k, p, conj(u[i]), v[i])
    end
    return k-1
end

function inner_product!(
    p::Vector{Term{T,Monomial{(x, y, z),3}}},
    u::vptype{T},
    v::vptype{T},
    cmat::Array{T,3}
    ) where T

    out = zero(eltype(cmat))
    iterm=dotpp!(p,u,v)
    @inbounds for ip=1:iterm
        c=coefficient(p[ip])
        exps=exponents(p[ip])
        out+=c*cmat[(exps .+ 1)...]
    end
    return out
end

#threaded and cached polynomials version of the force projection function
function projectforcet!(
    A::AbstractArray{T,2},
    cmat::Array{T,3},
    vs_i::Array{Array{P,1},1},
    vs_j::Array{Array{P,1},1},
    forcefun::Function, 
    args...; 
    n_cache = 10^6
    ) where {T, P <: Polynomial{T}}

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    @assert n_1 == size(A,1)
    @assert n_2 == size(A,2)
	# p = Progress(n_1*n_2)

    @sync for j = 1:n_2
        Threads.@spawn begin
           ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache) #cache terms array.
            f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
            for i = 1:n_1
                A[i,j] =  inner_product!(ptemp, vs_i[i], f, cmat)
            end
        end
    end
end
