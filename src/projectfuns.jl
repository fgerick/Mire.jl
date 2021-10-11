## low level functions that project the basis elements on a force function (forcefun)



# """
#     projectforce(N::Integer,cmat::Array{T,3},vs_i::Array{Array{P,1},1},vs_j::Array{Array{P,1},1},forcefun::Function,a::T,b::T,c::T, args...) where {T, P <: Polynomial{T}}

# Allocates new matrix `A` and fills elements by calling
# projectforce!(A,cmat,vs_i,vs_j,forcefun, args...)

# where `cmat[i,j,k]` contains the integrals of monomials xⁱyʲzᵏ.

# #Arguments:
# - `N`: maximum monomial degree
# - `vs_i`: basis vectors to project on
# - `vs_j`: basis vectors used for `forcefun`
# - `forcefun`: function of the force, e.g. coriolis
# - `args`: other arguments needed for `forcefun`
# """
# function projectforce(
#     cmat::Array{T,3},
#     vs_i::Array{Array{P,1},1},
#     vs_j::Array{Array{P,1},1},
#     forcefun::Function, 
#     args...;
#     kwargs...
#     ) where {T, P <: Polynomial{T}}

#     n_1 = length(vs_i)
#     n_2 = length(vs_j)

#     A = spzeros(T,n_1,n_2)

#     projectforce!(A, cmat, vs_i, vs_j, forcefun, args...; kwargs...)

#     return A
# end

# function projectforce(
#     cmat::Array{T,3},
#     vs::Array{Array{P,1},1},
#     forcefun::Function, 
#     args...;
#     kwargs...
#     ) where {T, P <: Polynomial{T}}

#     return projectforce(cmat,vs,vs,forcefun,args...; kwargs...)
# end

# """
#     projectforce!(A::AbstractArray{T, 2}, cmat::Array{T, 3}, vs::Array{Array{P, 1}, 1}, N::Integer, forcefun::Function, a::T, b::T, c::T, args...) where {T, P <: Polynomial{T}}

# DOCSTRING

# #Arguments:
# - `A`: pre-allocated array
# - `cmat`: pre-cached monomial integration values
# - `vs`: basis vectors
# - `N`: maximum monomial degree
# - `forcefun`: function of the force, e.g. coriolis
# - `args`: other arguments needed for `forcefun`
# """
# function projectforce!(
#     A::AbstractArray{T,2},
#     cmat::Array{T,3},
#     vs::Array{Array{P,1},1}, 
#     forcefun::Function, 
#     args...;
#     kwargs...
#     ) where {T, P <: Polynomial{T}}

#     return projectforce!(A, cmat, vs, vs, forcefun, args...; kwargs...)
# end


# function projectforce!(
#     A,
#     cmat, 
#     vs_i, 
#     vs_j, 
#     forcefun, 
#     args...; 
#     verbose=false,
#     thresh=10eps())

#     n_1 = length(vs_i)
#     n_2 = length(vs_j)
#     if verbose
#         p = Progress(n_1*n_2)
#     end

#     for j = 1:n_2
#         f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
#         for i = 1:n_1
#             aij = inner_product(vs_i[i], f, cmat)
#             if abs(aij) > thresh
#                 A[i,j] = aij
#             end
#             if verbose
#                 next!(p)
#             end
#         end
#     end
#     return nothing #vcat(itemps...), vcat(jtemps...), vcat(valtemps...)
# end

#vector type definition/abbreviation
# ptype{T} = Polynomial{T,Term{T,Monomial{(x, y, z),3}},Array{Term{T,Monomial{(x, y, z),3}},1}}
# vptype{T} = Vector{ptype{T}}


function projectforcet(vs_i, vs_j, cmat, forcefun, args...; kwargs...)
    nt = Threads.nthreads()
    itemps = [Int[] for i=1:nt]
    jtemps = [Int[] for i=1:nt]
    valtemps = [coefficienttype(vs_i[1][1])[] for i=1:nt]
     projectforcet!(0, 0, itemps, jtemps, valtemps, cmat, vs_i, vs_j, forcefun, args...; kwargs...)
    return sparse(vcat(itemps...), vcat(jtemps...), vcat(valtemps...), length(vs_i), length(vs_j))
end

function projectforcet!(
    i0,
    j0,
    itemps,
    jtemps,
    valtemps,
    cmat, 
    vs_i, 
    vs_j, 
    forcefun, 
    args...; 
    verbose=false,
    thresh=10eps())

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    if verbose
        p = Progress(n_1*n_2)
    end

    for j = 1:n_2
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        
        @sync for i = 1:n_1
            Threads.@spawn begin
                id = Threads.threadid()
                aij = inner_product(vs_i[i], f, cmat)
                if abs(aij) > thresh
                    push!(itemps[id], i+i0)
                    push!(jtemps[id], j+j0)
                    push!(valtemps[id], aij)
                    # A[i,j] = aij
                end
                if verbose
                    next!(p)
                end
            end
        end
    end
    return nothing #vcat(itemps...), vcat(jtemps...), vcat(valtemps...)
end
#
function projectforce(vs_i, vs_j, cmat, forcefun, args...; kwargs...)
    itemps = Int[]
    jtemps = Int[]
    valtemps = coefficienttype(vs_i[1][1])[]
     projectforce!(0, 0, itemps, jtemps, valtemps, cmat, vs_i, vs_j, forcefun, args...; kwargs...)
    return sparse(itemps,jtemps,valtemps, length(vs_i), length(vs_j))
end

function projectforce!(
    i0,
    j0,
    itemps,
    jtemps,
    valtemps,
    cmat, 
    vs_i, 
    vs_j, 
    forcefun, 
    args...; 
    verbose=false,
    thresh=10eps())

    n_1 = length(vs_i)
    n_2 = length(vs_j)
    if verbose
        p = Progress(n_1*n_2)
    end

    for j = 1:n_2
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        
        for i = 1:n_1
            aij = inner_product(vs_i[i], f, cmat)
            if abs(aij) > thresh
                push!(itemps, i+i0)
                push!(jtemps, j+j0)
                push!(valtemps, aij)
                # A[i,j] = aij
            end
            if verbose
                next!(p)
            end
        end
    end
    return nothing #vcat(itemps...), vcat(jtemps...), vcat(valtemps...)
end


#threaded and cached polynomials version of the force projection functioN
"""
    projectforcet!(args...)

Multithreaded projection 
"""
function projectforcet! end



function projectforcet_symmetric_neighbours!(
    i0,
    j0,
    itemps,
    jtemps,
    valtemps,
    cmat,
    vs_i,
    vs_j,
    forcefun::Function, 
    ls::Vector{Int},
    ms::Vector{Int},
    ispt::Vector{Bool}, 
    args...; 
    thresh=10eps(),
    )

    n_1 = length(vs_i)
    n_2 = length(vs_j)
 
    @sync for j=1:n_2
        # id = Threads.threadid()
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        for i = j:n_1
            if (ls[i]==ls[j]) && (ms[i]==ms[j]) && (ispt[i]==ispt[j])
                Threads.@spawn begin
                    id = Threads.threadid()
                    aij =  inner_product(vs_i[i], f, cmat)
                    if abs(aij) > thresh
                        push!(itemps[id],i+i0)
                        push!(jtemps[id],j+j0)
                        push!(valtemps[id],aij)
                        if i!=j
                            push!(itemps[id],j+j0)
                            push!(jtemps[id],i+i0)
                            push!(valtemps[id],aij)
                        end
                    end
                end
            end

        end
    end
end

function projectforce_symmetric_neighbours!(
    A,
    cmat,
    vs_i,
    vs_j,
    forcefun::Function, 
    ls::Vector{Int},
    ms::Vector{Int},
    ispt::Vector{Bool}, 
    args...; 
    thresh=10eps(),
    )

    n_1 = length(vs_i)
    n_2 = length(vs_j)
 
    for j=1:n_2
        f = forcefun(vs_j[j],args...) #calculate f(uⱼ)
        for i = j:n_1
            if (ls[i]==ls[j]) && (ms[i]==ms[j]) && (ispt[i]==ispt[j])
                aij =  inner_product(vs_i[i], f, cmat)
                if abs(aij) > thresh
                    A[i,j] = aij
                    if i!=j
                       A[j,i] = aij
                    end
                end
            end

        end
    end
end



# function projectlorentzqg!(i0, j0, itemps, jtemps, valtemps, cmat, vbasis, bbasis, b0, ls, ms, msu, lb0, mb0, ispoloidal; thresh=eps())

#     n_1 = length(vbasis)
#     n_2 = length(bbasis)
 
#     for j = 1:n_2
#         # id = Threads.threadid()
#         f = lorentz(bbasis[j],b0) #calculate ∇×Bⱼ×B₀ + ∇×B₀×Bⱼ
#         for i = 1:n_1
# 			if ispoloidal[j]
# 				if iseven(msu[i])
# 					lcondition = iseven(lb0) ? iseven(ls[j]) : isodd(ls[j])
# 				else
# 					lcondition = iseven(lb0) ? isodd(ls[j]) : iseven(ls[j])
# 				end
# 			else
# 				if iseven(msu[i])
# 					lcondition = iseven(lb0) ? isodd(ls[j]) : iseven(ls[j])
# 				else
# 					lcondition = iseven(lb0) ? iseven(ls[j]) : isodd(ls[j])
# 				end
# 			end	

# 			mcondition = ms[j] ∈ (msu[i]+mb0, msu[i]-mb0)
			 
#             if lcondition && mcondition
# 				aij =  inner_product(vbasis[i], f, cmat)
# 				if abs(aij) > thresh
# 					push!(itemps,i+i0)
# 					push!(jtemps,j+j0)
# 					push!(valtemps,aij)
# 				end
#             end

#         end
#     end
# end

function projectlorentzqgt!(i0, j0, itemps, jtemps, valtemps, cmat, vbasis, bbasis, b0, ls, ms, msu, lb0, mb0, ispoloidal, b0isp; thresh=eps())

    n_1 = length(vbasis)
    n_2 = length(bbasis)
 
    for j = 1:n_2
        
        f = lorentz(bbasis[j],b0) #calculate ∇×Bⱼ×B₀ + ∇×B₀×Bⱼ
        for i = 1:n_1
			if (b0isp ? ispoloidal[j] : !ispoloidal[j])
				if iseven(msu[i])
					lcondition = iseven(lb0) ? iseven(ls[j]) : isodd(ls[j])
				else
					lcondition = iseven(lb0) ? isodd(ls[j]) : iseven(ls[j])
				end
			else
				if iseven(msu[i])
					lcondition = iseven(lb0) ? isodd(ls[j]) : iseven(ls[j])
				else
					lcondition = iseven(lb0) ? iseven(ls[j]) : isodd(ls[j])
				end
			end	

			mcondition = abs(ms[j]) ∈ (abs(msu[i]+mb0), abs(msu[i]-mb0))
			 
            if lcondition && mcondition
                Threads.@spawn begin
                    id = Threads.threadid()
                    aij =  inner_product(vbasis[i], f, cmat)
                    if abs(aij) > thresh
                        push!(itemps[id],i+i0)
                        push!(jtemps[id],j+j0)
                        push!(valtemps[id],aij)
                    end
                end
            end

        end
    end
end

# function projectlorentzqg(vs_i, vs_j, cmat, args...; kwargs...)
# 	# nt = Threads.nthreads()
# 	itemps = Int[]
# 	jtemps = Int[]
# 	valtemps = coefficienttype(vs_i[1][1])[]
# 		projectlorentzqg!(0, 0, itemps, jtemps, valtemps, cmat, vs_i, vs_j, args...; kwargs...)
# 	return sparse(vcat(itemps...), vcat(jtemps...), vcat(valtemps...), length(vs_i), length(vs_j))
# end

# function projectinductionqg!(i0, j0, itemps, jtemps, valtemps, cmat, bbasis, vbasis, b0, ls, ms, msu, lb0, mb0, ispoloidal; thresh=eps())

#     n_1 = length(bbasis)
#     n_2 = length(vbasis)
 
#     for j = 1:n_2
#         # id = Threads.threadid()
#         f = Mire.advection(vbasis[j],b0) #calculate ∇×(uⱼ×b0)
#         for i = 1:n_1
# 			if ispoloidal[i]
# 				if iseven(msu[j])
# 					lcondition = iseven(lb0) ? iseven(ls[i]) : isodd(ls[i])
# 				else
# 					lcondition = iseven(lb0) ? isodd(ls[i]) : iseven(ls[i])
# 				end
# 			else
# 				if iseven(msu[j])
# 					lcondition = iseven(lb0) ? isodd(ls[i]) : iseven(ls[i])
# 				else
# 					lcondition = iseven(lb0) ? iseven(ls[i]) : isodd(ls[i])
# 				end
# 			end	

# 			mcondition = ms[i] ∈ (msu[j]+mb0, msu[j]-mb0)
			 
#             if lcondition && mcondition
# 				aij =  inner_product(bbasis[i], f, cmat)
# 				if abs(aij) > thresh
# 					push!(itemps,i+i0)
# 					push!(jtemps,j+j0)
# 					push!(valtemps,aij)
# 				end
#             end

#         end
#     end
# end

function projectinductionqgt!(i0, j0, itemps, jtemps, valtemps, cmat, bbasis, vbasis, b0, ls, ms, msu, lb0, mb0, ispoloidal, b0isp; thresh=eps())

    n_1 = length(bbasis)
    n_2 = length(vbasis)
 
    @sync for j = 1:n_2
        f = Mire.advection(vbasis[j],b0) #calculate ∇×(uⱼ×b0)
        for i = 1:n_1
			if (b0isp ? ispoloidal[i] : !ispoloidal[i])
				if iseven(msu[j])
					lcondition = iseven(lb0) ? iseven(ls[i]) : isodd(ls[i])
				else
					lcondition = iseven(lb0) ? isodd(ls[i]) : iseven(ls[i])
				end
			else
				if iseven(msu[j])
					lcondition = iseven(lb0) ? isodd(ls[i]) : iseven(ls[i])
				else
					lcondition = iseven(lb0) ? iseven(ls[i]) : isodd(ls[i])
				end
			end	

			mcondition = abs(ms[i]) ∈ (abs(msu[j]+mb0), abs(msu[j]-mb0))
			 
            if lcondition && mcondition
                Threads.@spawn begin
                    id = Threads.threadid()
                    aij =  inner_product(bbasis[i], f, cmat)
                    if abs(aij) > thresh
                        push!(itemps[id],i+i0)
                        push!(jtemps[id],j+j0)
                        push!(valtemps[id],aij)
                    end
                end
            end

        end
    end
end

function projectcoriolisqgt!(i0, j0, itemps, jtemps, valtemps, cmat, vbasis, Ω; thresh=eps())

    n_1 = length(vbasis)
 
    Threads.@threads for i = 1:n_1
        id = Threads.threadid()
        f = Mire.coriolis(vbasis[i],Ω) #calculate 2uⱼ×Ω
        aij =  inner_product(vbasis[i], f, cmat)
        if abs(aij) > thresh
            push!(itemps[id],i+i0)
            push!(jtemps[id],i+j0)
            push!(valtemps[id],aij)
        end
    end
    return nothing
end
# function projectinductionqg(vs_i, vs_j, cmat, args...; kwargs...)
# 	# nt = Threads.nthreads()
# 	itemps = Int[]
# 	jtemps = Int[]
# 	valtemps = coefficienttype(vs_i[1][1])[]
# 		projectinductionqg!(0, 0, itemps, jtemps, valtemps, cmat, vs_i, vs_j, args...; kwargs...)
# 	return sparse(vcat(itemps...), vcat(jtemps...), vcat(valtemps...), length(vs_i), length(vs_j))
# end
