
#pseudo vacuum BC (Vidal & Cebron, 2021)

struct PVMFBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end

function normal(V::Sphere{T}) where T
	return [x,y,z]
end

function normal(V::Ellipsoid{T}) where T
	return [x/V.a^2,y/V.b^2,z/V.c^2]
end

function poincare(p,V::Ellipsoid{T}) where T
	return p(x=>x/V.a,y=>y/V.b,z=>z/V.c)
end

poincare(p, V::Sphere{T}) where T = p
δ(i,j,k) = i==j==k

function ℬ(Φ,V::Sphere{T}) where T
	ΔΦ = Δ(Φ)
	ℬ = zero(Φ)
	for (m,c) in zip(monomials(ΔΦ),coefficients(ΔΦ))
		i,j,k = exponents(m)
		ℬ += -c*m/(i+j+k+3)
	end
	return ℬ
end

function ℬ(Φ,V::Ellipsoid{T}) where T
	ΔΦ = Δ(Φ)
	ℬ = zero(Φ)
	for (m,c) in zip(monomials(ΔΦ),coefficients(ΔΦ))
		i,j,k = exponents(m)
		ℬ += -c*m/((i+1)/V.a^2+(j+1)/V.b^2+(k+1)/V.c^2)
	end
	return ℬ
end

function basisvector_A(::Type{PVMFBasis}, V::Volume{T}, l, m, p; kwargs...) where T
	n = normal(V)
	rlm_ = poincare(rlm(l,m,x,y,z; kwargs...),V)
	𝒜 = F(V)*(1-F(V))^p*rlm_
	return ∇ × (𝒜*n)
end

function basisvector_B(::Type{PVMFBasis}, V::Volume{T}, l, m, p; kwargs...) where T
	n = normal(V)
	rlm_ = poincare(rlm(l,m,x,y,z; kwargs...),V)
	Φ = F(V)*(1-F(V))^p*rlm_
	return ℬ(Φ,V)*n + ∇(Φ)
end

function basisvectors(::Type{PVMFBasis}, N::Int, V::Volume{T}; m=-N:N, kwargs...) where T

	b_A = [basisvector_A(PVMFBasis, V, l,m,p; kwargs...) for l=1:(N-2) for m in m for p=(0:(N-l)÷2-1) if abs(m)<=l]
	b_B = [basisvector_B(PVMFBasis, V, l,m,p; kwargs...) for l=1:(N-1) for m in m for p=(0:(N+1-l)÷2-1) if abs(m)<=l]
	return vcat(b_A, b_B)
end

PVMFBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} =
                PVMFBasis(N, m, V, basisvectors(PVMFBasis, N, V; m, kwargs...),false)