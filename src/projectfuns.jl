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
    args...;
    kwargs...
    ) where {T, P <: Polynomial{T}}

    return projectforce!(A, cmat, vs, vs, forcefun, args...; kwargs...)
end

# function projectforce!(
#     A::AbstractArray{T,2},
#     cmat::Array{T,3},
#     vs_i::Array{Array{P,1},1},
#     vs_j::Array{Array{P,1},1},
#     forcefun::Function,
#     args...
#     ) where {T, P <: Polynomial{T}}

#     n_1 = length(vs_i)
#     n_2 = length(vs_j)
#     @assert n_1 == size(A,1)
#     @assert n_2 == size(A,2)

#     @inbounds for j=1:n_2
#         f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
#         for i=1:n_1
#             A[i,j] = inner_product_real(cmat,vs_i[i],f)
#         end
#     end

#     return nothing
# end


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
    args...;
    kwargs...
    ) where {T, P <: Polynomial{T}}

    n_1 = length(vs_i)
    n_2 = length(vs_j)

    A = spzeros(T,n_1,n_2)

    projectforce!(A, cmat, vs_i, vs_j, forcefun, args...; kwargs...)

    return A
end

function projectforce(
    cmat::Array{T,3},
    vs::Array{Array{P,1},1},
    forcefun::Function, 
    args...;
    kwargs...
    ) where {T, P <: Polynomial{T}}

    return projectforce(cmat,vs,vs,forcefun,args...; kwargs...)
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
            p[iter] = conj(coefficient(ti))*monomial(ti)*tj
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
        k = _mul_to_terms!(k, p, u[i], v[i])
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

#threaded and cached polynomials version of the force projection functioN
"""
    projectforcet!(A::Array{T,2}, cmat::Array{T,3}, vs_i::Array{Array{P,1},1}, vs_j::Array{Array{P,1},1}, forcefun::Function, args...; n_cache = 10^6) where {T, P <: Polynomial{T}}

Multithreaded projection using pre-cached polynomial terms array. `n_cache` needs to be sufficiently large for a given
polynomial degree of the basis vectors (trial and error for now!). Matrix `A` must be dense due to thread safety.
"""
# function projectforcet!(
#     A::AbstractArray{T,2},
#     cmat::Array{T,3},
#     vs_i::Array{Array{P,1},1},
#     vs_j::Array{Array{P,1},1},
#     forcefun::Function, 
#     args...; 
#     n_cache = 10^6
#     ) where {T, P <: Polynomial{T}}

#     n_1 = length(vs_i)
#     n_2 = length(vs_j)
#     @assert n_1 == size(A,1)
#     @assert n_2 == size(A,2)
# 	# p = Progress(n_1*n_2)

#     @sync for j = 1:n_2
#         Threads.@spawn begin
#            ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache) #cache terms array.
#             f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
#             for i = 1:n_1
#                 A[i,j] =  inner_product!(ptemp, vs_i[i], f, cmat)
#             end
#         end
#     end
# end

function projectforcet!(
    A::AbstractMatrix{T},
    cmat::Array{T,3},
    vs_i::Union{Vector{vptype{T}},Vector{vptype{T2}}},
    vs_j::Union{Vector{vptype{T}},Vector{vptype{T2}}},
    forcefun::Function, 
    args...; 
    n_cache = 10^6,
    verbose = false
    ) where {T, T2}

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    @assert n_1 == size(A,1)
    @assert n_2 == size(A,2)
    if verbose
        p = Progress(n_1*n_2)
    end
    nt = Threads.nthreads()
    ptemps = [zeros(Term{T,Monomial{(x, y, z),3}}, n_cache) for i=1:nt]
    # ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache)
    # Threads.@threads 
    Threads.@threads for j = 1:n_2
        #    ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache) #cache terms array.
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        for i = 1:n_1
            A[i,j] =  inner_product!(ptemps[Threads.threadid()], vs_i[i], f, cmat)
            if verbose
                next!(p)
            end
        end
    end
end


function projectforce!(
    A::AbstractMatrix{T},
    cmat::Array{T,3},
    vs_i::Union{Vector{vptype{T}},Vector{vptype{T2}}},
    vs_j::Union{Vector{vptype{T}},Vector{vptype{T2}}},
    forcefun::Function, 
    args...; 
    n_cache = 10^6,
    verbose = false
    ) where {T, T2}

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    @assert n_1 == size(A,1)
    @assert n_2 == size(A,2)
    nt = Threads.nthreads()
    ptemp = zeros(Term{T,Monomial{(x, y, z),3}}, n_cache)
    if verbose
        p = Progress(n_1*n_2)
    end
    @inbounds for j = 1:n_2
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        for i = 1:n_1
            A[i,j] =  inner_product!(ptemp, vs_i[i], f, cmat)
            if verbose
                next!(p)
            end
        end
    end
end

function projectforcet_symmetric_neighbours!(
    A::AbstractMatrix{T},
    cmat::Array{T,3},
    vs_i::Union{Vector{vptype{T}},Vector{vptype{T2}}},
    vs_j::Union{Vector{vptype{T}},Vector{vptype{T2}}},
    forcefun::Function, 
    ls::Vector{Int},
    ms::Vector{Int},
    ispt::Vector{Bool}, 
    args...; 
    n_cache = 2*10^6,
    verbose = false
    ) where {T, T2}

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    @assert n_1 == size(A,1)
    @assert n_2 == size(A,2)
    # ptemp = zeros(Term{T,Monomial{(x, y, z),3}},n_cache);
    nt = Threads.nthreads()
    ptemps = [zeros(Term{T,Monomial{(x, y, z),3}}, n_cache) for i=1:nt]

    Threads.@threads for j=1:n_2
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        for i = j:n_1
            if (ls[i]==ls[j]) && (ms[i]==ms[j]) && (ispt[i]==ispt[j])
                A_ij =  inner_product!(ptemps[Threads.threadid()], vs_i[i], f, cmat)
                A[i,j] = A_ij
                if i!=j
                    A[j,i] = A_ij
                end
            end

        end
    end
end
