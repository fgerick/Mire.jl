
## Basis following Ivers (2017) 

L(V::Ellipsoid) = [1/V.a 0 0
			0 1/V.b 0
			0 0 1/V.c]

L(V::Sphere) = I

radius(V::Sphere) = V.r₀
radius(V::Ellipsoid) = 1 #todo for the moment

function s(l,m,n,V::Volume; kwargs...)
	∇̇ = Del(V)
	Plm = poincare((radius(V)^2-r²)*r²^n*rlm(l,m,x,y,z; kwargs...),V)
	𝐫̇ = L(V)*𝐫
	return inv(L(V))*(∇̇×(∇̇×(Plm*𝐫̇)))
end

function t(l,m,n,V::Volume; kwargs...)
	∇̇ = Del(V)
	Tlm = poincare(r²^n*rlm(l,m,x,y,z; kwargs...),V)
	𝐫̇ = L(V)*𝐫
	return inv(L(V))*(∇̇×(Tlm*𝐫̇))
end

struct IversBasisC{T<:Number,Vol<:Volume{T}} <: Mire.VectorBasis{T,Vol}
	N::Int
	m
	V::Vol
	el::Vector{Mire.vptype{Complex{T}}}
	orthonorm::Bool
end

function basisvectors(::Type{IversBasisC}, N::Int, V::Volume{T}; m=-N:N, norm=Mire.Schmidt{T}) where T
	lmn_p = [(l,m,n) for l = 1:(N-1) for m in m for n = 0:((N+1-l)÷2-1) if abs(m)<=l]
	lmn_t = [(l,m,n) for l = 1:N for m in m for n = 0:((N-l)÷2) if abs(m)<=l]
	u_p = [s(l,m,n,V; real=false, norm) for (l,m,n) in lmn_p]
	u_t = [t(l,m,n,V; real=false, norm) for (l,m,n) in lmn_t]
	
	return vcat(u_p,u_t)
end
	

function LMN(B::IversBasisC)
	N = B.N
	lmn_p = [(l,m,n) for l = 1:(N-1) for m = 0:l for n = 0:((N+1-l)÷2-1)]
	lmn_t = [(l,m,n) for l = 1:N for m = 0:l for n = 0:((N-l)÷2)]
	return lmn_p,lmn_t
end
	
struct IversBasis{T<:Number,Vol<:Volume{T}} <: Mire.VectorBasis{T,Vol}
	N::Int
	m
	V::Vol
	el::Vector{Mire.vptype{T}}
	orthonorm::Bool
end

function basisvectors(::Type{IversBasis}, N::Int, V::Volume{T}; norm=Mire.Schmidt{T}, m=-N:N) where T
	lmn_p = [(l,m,n) for l = 1:(N-1) for m in m for n = 0:((N+1-l)÷2-1) if abs(m)<=l]
	lmn_t = [(l,m,n) for l = 1:N for m in m for n = 0:((N-l)÷2) if abs(m)<=l]
	u_p = [s(l,m,n,V; real=true, norm) for (l,m,n) in lmn_p]
	u_t = [t(l,m,n,V; real=true, norm) for (l,m,n) in lmn_t]
	
	return vcat(u_p,u_t)
end
	

function LMN(B::IversBasis)
	N = B.N
	m = B.m
	lmn_p = [(l,m,n) for l = 1:(N-1) for m in m for n = 0:((N+1-l)÷2-1) if abs(m)<=l]
	lmn_t = [(l,m,n) for l = 1:N for m in m for n = 0:((N-l)÷2) if abs(m)<=l]
	return lmn_p,lmn_t
end

IversBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} =
                IversBasis(N, m, V, basisvectors(IversBasis, N, V; m, kwargs...),false)

IversBasisC(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} =
                IversBasisC(N, m, V, basisvectors(IversBasisC, N, V; m, kwargs...),false)