
#MAGNETIC FIELD BASES:


"""
	ConductingMFBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}

Basis of 3-D magnetic field, so that \$\\mathbf{B}\\cdot\\mathbf{n} = 0\$ at \$\\partial\\mathcal{V}\$
and \$\\nabla\\cdot\\mathbf{B} = 0\$ after Lebovitz (1989). It is exactly the same as `LebovitzBasis`.
"""
const ConductingMFBasis = LebovitzBasis

"""
InsulatingMFBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}

Basis of insulating magnetic fields following Gerick et al. (2021). For now, only for `Vol<:Sphere{T}`, i.e. in a spherical domain. 
"""
struct InsulatingMFBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end

#coefficients of interior and exterior potential field correction
vi(l::Int, n::Int) = -(l + 1)//(4l*n + 4l + 2n + 2)
ve(l::Int, n::Int) = l//((2*l + 1)*(2l + 2n + 3))

function bpol_g21(l,m,n; kwargs...)
	Rₗᵐ = rlm(l,m,x,y,z; kwargs...)
	Pₗₘₙ = -r²^(n+1)*Rₗᵐ / (2(n+1)*(2l+2n+3))
	return ∇×(∇×(Pₗₘₙ*𝐫)) - vi(l,n)*∇(Rₗᵐ)
end

function btor_g21(l,m,n; kwargs...)
	Rₗᵐ = rlm(l,m,x,y,z; kwargs...)
	Tₗₘₙ = (1-r²)*r²^n*Rₗᵐ
	return ∇×(Tₗₘₙ*𝐫)
end

function LMN(P::InsulatingMFBasis)
    N = P.N
	m = P.m

	lmn_p = [(l,m,n) for l = 1:(N-1) for m in m for n = 0:(N-l+1)÷2-1 if abs(m)<=l]
	lmn_t = [(l,m,n) for l = 1:(N-2) for m in m for n = 0:(N-l)÷2-1 if abs(m)<=l]
	
	return lmn_p,lmn_t
end

function basisvectors(::Type{InsulatingMFBasis}, N::Int, V::Sphere{T}; m=-N:N, norm=Schmidt{T}) where T
    if typeof(V) != Sphere{T}
        return throw(ArgumentError("Insulating magnetic field basis is only implemented in the sphere!"))
    end

	lmn_p = [(l,m,n) for l = 1:(N-1) for m in m for n = 0:(N-l+1)÷2-1 if abs(m)<=l]
	lmn_t = [(l,m,n) for l = 1:(N-2) for m in m for n = 0:(N-l)÷2-1 if abs(m)<=l]
	
	bpols = [bpol_g21(lmn...; real=true, norm) for lmn in lmn_p]
	btors = [btor_g21(lmn...; real=true, norm) for lmn in lmn_t]

	return vcat(bpols,btors)
end

InsulatingMFBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = InsulatingMFBasis(N, m, V, basisvectors(InsulatingMFBasis, N, V; m, kwargs...),false)

"""
InsulatingMFCBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}

Complex basis of insulating magnetic fields following Gerick et al. (2021). For now, only for `Vol<:Sphere{T}`, i.e. in a spherical domain. 
"""
struct InsulatingMFCBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{Complex{T}}}
    orthonorm::Bool
end


function basisvectors(::Type{InsulatingMFCBasis}, N::Int, V::Sphere{T}; m=-N:N, norm=Schmidt{T}) where T
    if typeof(V) != Sphere{T}
        return throw(ArgumentError("Insulating magnetic field basis is only implemented in the sphere!"))
    end

	lmn_p = [(l,m,n) for l = 1:(N-1) for m = m for n = 0:(N-l+1)÷2-1 if abs(m)<=l]
	lmn_t = [(l,m,n) for l = 1:(N-2) for m = m for n = 0:(N-l)÷2-1 if abs(m)<=l]
	
	bpols = [bpol_g21(lmn...; real=false, norm) for lmn in lmn_p]
	btors = [btor_g21(lmn...; real=false, norm) for lmn in lmn_t]

	return vcat(bpols,btors)
end

function LMN(P::InsulatingMFCBasis)
    N = P.N

	lmn_p = [(l,m,n) for l = 1:(N-1) for m = 0:l for n = 0:(N-l+1)÷2-1]
	lmn_t = [(l,m,n) for l = 1:(N-2) for m = 0:l for n = 0:(N-l)÷2-1]
	
	return lmn_p,lmn_t
end


InsulatingMFCBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = InsulatingMFCBasis(N, m, V, basisvectors(InsulatingMFCBasis, N, V; m, kwargs...),false)

struct InsMFONCBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{Complex{T}}}
    orthonorm::Bool
end


h(::Type{T}, l::Integer, n::Integer) where T = (one(T)-r2)*jacobi(2r2-1, n-1, T(2), T(l + 1//2))

function k(::Type{T}, l::Integer, n::Integer) where T
	c₀ = T(-2n^2*(l + 1) - n*(l + 1)*(2l - 1) - l*(2l + 1))
	c₁ = T(2*(l + 1)*n^2 + (2l + 3)*(l + 1)*n + (2l + 1)^2)
	c₂ = T(4n*l + l*(2l+1))
	P₀ = jacobi(2r2-1, n, zero(T), T(l + 1//2))
	P₁ = jacobi(2r2-1, n-1, zero(T), T(l + 1//2))
	return c₀*P₀ + c₁*P₁ + c₂
end

btor(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (h(T,l,n)*rlm(l,m,x,y,z; norm=Schmidt{T}, kwargs...)*𝐫)
bpol(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (∇ × (k(T,l,n)*rlm(l,m,x,y,z; norm=Schmidt{T}, kwargs...)*𝐫))


function basisvectors(::Type{InsMFONCBasis}, N::Int, V::Sphere{T}; m=-N:N, kwargs...) where T
	N-=1

	lmnsp = [(l,m,n)  for l in 1:N for m in m for n in 1:(N-l+2)÷2 if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:(N-1) for m in m for n in 1:((N+1-l)÷2) if abs(m)<=l]

	BP = map(lmn->bpol(T,lmn...; real=false, kwargs...),lmnsp)
	BT = map(lmn->btor(T,lmn...; real=false, kwargs...),lmnst)
    return vcat(BP,BT)
	# return BP, BT, ls, ms, ns, lstor,mstor,nstor
end



function LMN(P::InsMFONCBasis{T,V}) where {T,V}
	# r2 = x^2 + y^2 + z^2
    N = P.N-1

	lmnsp = [(l,m,n)  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]

	return lmnsp, lmnst
end


InsMFONCBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = InsMFONCBasis(N, m, V, basisvectors(InsMFONCBasis, N, V; m, kwargs...),false)

struct InsMFONBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end


function basisvectors(::Type{InsMFONBasis}, N::Int, V::Sphere{T}; m=-N:N, kwargs...) where T
	# N-=1

	lmnsp = [(l,m,n)  for l in 1:(N-1) for m in m for n in 1:(N-l+1)÷2 if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:(N-2) for m in m for n in 1:((N-l)÷2) if abs(m)<=l]

	BP = map(lmn->bpol(T,lmn...; kwargs...),lmnsp)
	BT = map(lmn->btor(T,lmn...; kwargs...),lmnst)
    return vcat(BP,BT)
end


function LMN(P::InsMFONBasis{T,V}) where {T,V}
    N = P.N
	m = P.m
	lmnsp = [(l,m,n)  for l in 1:(N-1) for m in m for n in 1:(N-l+1)÷2 if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:(N-2) for m in m for n in 1:((N-l)÷2) if abs(m)<=l]
	
	return lmnsp, lmnst
end



InsMFONBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = InsMFONBasis(N, m, V, basisvectors(InsMFONBasis, N, V; m, kwargs...),false)


function t_mf_lapl(T,l,n) 
    fac = 1/sqrt(l*(1 + l)*(1/T(-1 + 2*l + 4*n) + 1/T(3 + 2*l + 4*n)))
    return fac * (jacobi(2r2-1,n,T(0),T(l+1//2)) - jacobi(2r2-1,n-1,T(0),T(l+1//2))) 
end


function s_mf_lapl(T,l, n)
    fac = 1/(sqrt(2l*T(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
    return fac * ( (2*l + 4*n - 3) * jacobi(2r2-1,n,T(0),T(l+1//2)) - 2*(2*l + 4*n - 1)*jacobi(2r2-1,n-1,T(0),T(l+1//2)) +(2*l + 4*n + 1)*jacobi(2r2-1,n-2,T(0),T(l+1//2)))
end

btor_lapl(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (t_mf_lapl(T,l,n)*rlm(l,m,x,y,z; norm=Full{T}, kwargs...)*𝐫)
bpol_lapl(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (∇ × (s_mf_lapl(T,l,n)*rlm(l,m,x,y,z; norm=Full{T}, kwargs...)*𝐫))



struct InsMFLaplBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end


function basisvectors(::Type{InsMFLaplBasis}, N::Int, V::Sphere{T}; m=-N:N, kwargs...) where T
	lmnsp, lmnst = LMN(InsMFLaplBasis, N, m)
	BP = map(lmn->bpol_lapl(T,lmn...; kwargs...),lmnsp)
	BT = map(lmn->btor_lapl(T,lmn...; kwargs...),lmnst)
    return vcat(BP,BT)
end


function LMN(P::InsMFLaplBasis{T,V}) where {T,V}
	N = P.N
	m = P.m
	lmnsp, lmnst = LMN(InsMFLaplBasis, N, m)
	return lmnsp, lmnst
end

function LMN(::Type{InsMFLaplBasis}, N, m=-N:N) 
	lmnsp = [(l,m,n)  for l in 1:N for m in m for n in 1:(N-l+1)÷2 if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:N for m in m for n in 1:((N-l)÷2) if abs(m)<=l]
	return lmnsp, lmnst
end


InsMFLaplBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = InsMFLaplBasis(N, m, V, basisvectors(InsMFLaplBasis, N, V; m, kwargs...),false)