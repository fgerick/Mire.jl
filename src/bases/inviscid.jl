struct InviscidBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end

#inviscid 3D basis, sphere.
function t_in(T,l,n) 
    fac = sqrt(T(3+2l+4n))/sqrt(T(l*(l+1)))
    return jacobi(2r2-1,n,T(0),T(l+1//2))*fac
end

function s_in(T,l,n) 
    fac = sqrt(T(5+2l+4n))/sqrt(T(4l*(l+1)*(n+1)^2))
    return (1-r2)*jacobi(2r2-1,n,T(1),T(l+1//2))*fac
end

utor_in(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (t_in(real(T),l,n)*rlm(l,m,x,y,z; norm=Full{T}, real=false)*𝐫)
upol_in(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (∇ × (s_in(real(T),l,n)*rlm(l,m,x,y,z; norm=Full{T}, real=false)*𝐫))



function basisvectors(::Type{InviscidBasis}, N::Int, V::Sphere{T}; m=-N:N, kwargs...) where T
	lmnsp, lmnst = LMN(InviscidBasis, N, m)

	BP = map(lmn->upol_in(T,lmn...; kwargs...),lmnsp)
	BT = map(lmn->utor_in(T,lmn...; kwargs...),lmnst)
    return vcat(BP,BT)
end


function LMN(P::InviscidBasis{T,V}) where {T,V}
	N = P.N
	m = P.m
	lmnsp, lmnst = LMN(InviscidBasis, N, m)

	return lmnsp, lmnst
end

function LMN(::Type{InviscidBasis}, N, m=-N:N)
	lmnsp = [(l,m,n)  for l in 1:N for m in m for n in 0:((N-l+1)÷2-1) if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:N for m in m for n in 0:((N-l)÷2) if abs(m)<=l]
	return lmnsp, lmnst
end

InviscidBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} =
                InviscidBasis(N, m, V, basisvectors(InviscidBasis, N, V; m, kwargs...),true)