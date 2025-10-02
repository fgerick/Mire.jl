
ptype{T} = Polynomial{T,Term{T,Monomial{(x, y, z),3}},Array{Term{T,Monomial{(x, y, z),3}},1}}
vptype{T} = Vector{ ptype{T}}

"""
    VectorBasis{T,V}

Defines `abstract type` of a vector basis to for a `T,V` so that `V=Volume{T}`.
"""
abstract type VectorBasis{T,V} end




"""
    LebovitzBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}

Basis of 3-D vector field, so that \$\\mathbf{u}\\cdot\\mathbf{n} = 0\$ at \$\\partial\\mathcal{V}\$
and \$\\nabla\\cdot\\mathbf{u} = 0\$ after Lebovitz (1989).
"""
struct LebovitzBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end

### Lebovitz (1989) basis

# Monomials
Π(n::Integer, m::Integer, l::Integer) = x^n * y^m * z^l


# Lebovitz 1989, eq. (39b)
F(V::Ellipsoid{T}) where {T} = (1 - x^2 / V.a^2 - y^2 / V.b^2 - z^2 / V.c^2)
F(V::Sphere{T}) where {T} = (one(T) - x^2 - y^2 - z^2)

# Basis vectors, Lebovitz 1989, eq. (41-42)
v1(n::Int, m::Int, l::Int, V::Volume{T}) where {T} = ∇(Π(n, m, l) * F(V)) × ex
v2(n::Int, m::Int, l::Int, V::Volume{T}) where {T} = ∇(Π(n, m, l) * F(V)) × ey
v3(n::Int, m::Int, l::Int, V::Volume{T}) where {T} = ∇(Π(n, m, l) * F(V)) × ez


#Computes all combos i,j,k satisfying i+j+k < N. Returns two `Vector{Vector{Int}}` with
#the first one all `[i,j,k]` with `k=0`, and the second one with `1<=k<N`.
function combos(N::Int)
    [[i, j, k] for i = 0:N for j = 0:N for k = 0 if (i + j + k < N)],
    [[i, j, k] for i = 0:N for j = 0:N for k = 1:N if (i + j + k < N)]
end


# Compute all velocity basis vectors for a given maximal degree `N` (Lebovitz 1989, eq. 41-42)
function basisvectors(::Type{LebovitzBasis}, N::Int, V::Volume{T}) where {T}
    gp, hp = combos(N)
    v_1 = [v1(h..., V) for h in vcat(gp, hp)]
    v_2 = [v2(h..., V) for h in vcat(gp, hp)]
    v_3 = [v3(g..., V) for g in gp]

    return vcat(v_1, v_2, v_3)
end





"""
    QGBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}

Basis of QG vector field (Gerick et al., 2020), so that \$\\mathbf{u}=\\nabla(h^3x^ny^m)\\times\\nabla(z/h)\$.
"""
struct QGBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end

qg_combos(N::Integer) = [[i, j] for i = 0:N for j = 0:N if (i + j < N)]
qg_combos_sorted(N::Integer) = (N==0) ? [[0,0]] : vcat(qg_combos2(N-1),[[i,j] for i=0:N for j=0:N if (N-1<i+j<=N)])

#Generate basis vector \$\\mathbf{u}=\\nabla(h^3x^ny^m)\\times\\nabla(z/h)\$
function uqg(n::Integer, m::Integer, V::Volume{T}) where {T}

    if typeof(V) <: Sphere
        a, b, c = one(T), one(T), one(T)
    elseif typeof(V) <: Ellipsoid
        a, b, c = V.a, V.b, V.c
    else
        error("QG velocity not implemented for this Volume")
    end
    Π = one(T) * x^n * y^m * z^0
    h2 = c^2 * (1 - x^2*1/a^2 - y^2*1/b^2)
    ez = Mire.ez #[zero(Π),zero(Π),one(Π)]
    h∇h = [-c^2/a^2 * x, -c^2/b^2 * y, 0]
    return h2 * ∇(Π) × ez + 3 * Π * h∇h × ez - z * ∇(Π) × h∇h
end


function basisvectors(::Type{QGBasis}, N::Int, V::Volume{T}) where {T}
    cs = qg_combos(N)
    return [uqg(ci..., V) for ci in cs]
end


"""
    QGIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}

Basis of complex QG inertial modes (Maffei et al., 2017). This basis is orthonormal.
"""
struct QGIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{Complex{T}}}
    orthonorm::Bool
end

CSH=CartesianSphericalHarmonics

function jacobi(x, n::Integer, a::R, b::R) where R <: Number
    ox = one(x*a*n)
    zx = zero(x*a*n)
    if n==0
        return ox
    elseif n==1
        return ox*1//2 * (a - b + (a + b + 2)*x)
    end

    p0 = ox
    p1 = ox*1//2 * (a - b + (a + b + 2)*x)
    p2 = zx;

    for i = 1:(n-1)
		a1 = 2*(i+1)*(i+a+b+1)*(2*i+a+b)
		a2 = (2*i+a+b+1)*(a*a-b*b);
		a3 = (2*i+a+b)*(2*i+a+b+1)*(2*i+a+b+2);
		a4 = 2*(i+a)*(i+b)*(2*i+a+b+2);
		p2 = ox/a1 * ( (a2 + a3*x)*p1 - a4*p0);

        p0 = p1
        p1 = p2
    end

    return p2
end

const r2 = x^2 + y^2 + z^2

function convert_polynom(S::Type{T},p) where T
	c = coefficients(p)
	m = monomials(p)
	return sum(S.(c).*m)
end


function prefac(n::BigInt,m::BigInt)
	m = abs(m)
	ω = -1/(n*(2n+2m+1)+m/2+m^2/6)
	fac = gamma(n+3//2)*gamma(n+m)/(gamma(n)*(2n + 1//2 + m)*gamma(n + 3//2 + m))
	return -4big(pi)/ω*fac
end

function basiselementc(n::Integer,m::Integer,V::Sphere{T}) where T
	h² = 1-x^2-y^2
	ez =   [0,0,1]
	h∇h =   [-x,-y,0]
	s² = x^2 + y^2
	J = jacobi(2s²-1,big(n)-1,big(3)//2, big(abs(m))//1 )

	ψp =   J * ((-1)^(m<0)*im*CSH.sinsinpoly(big(abs(m)),x,y) + CSH.cossinpoly(big(abs(m)),x,y))
	# if m<0
	# 	ψp = conj(ψp)
	# end
	ψp = convert_polynom(complex(T),ψp/sqrt(prefac(big(n),big(m))))

	return h²*∇(ψp)×ez + 3*ψp*h∇h×ez - z*∇(ψp)×h∇h
end

function basisvectors(::Type{QGIMBasis}, N::Int, V::Volume{T}; m=-(N-1):(N-1)) where T
    return [basiselementc(n,m,V) for m in m for n=1:(N-abs(m)+1)÷2 ]
end


function NM(b::QGIMBasis; m = -(b.N-1):(b.N-1))
	N = b.N
	ms = [m for m in m for n=1:(N-abs(m)+1)÷2 ]
	ns = [n for m in m for n=1:(N-abs(m)+1)÷2 ]
	return ns,ms
end



"""
    QGRIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}

Real basis similar to the QG inertial modes (Maffei et al., 2017). They are not the solutions to the QG inertial mode equation and they are not orthogonal.
"""
struct QGRIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end


function basisvectors(::Type{QGRIMBasis}, N::Int, V::Volume{T}; m=-(N-1):(N-1)) where T
    return [basiselement(n,m,V) for m in m for n=1:(N-abs(m)-1)÷2]
end

function basiselement(n::Integer,m::Integer,V::Sphere{T}) where T
	h² = 1-x^2-y^2
	ez =   [0,0,1]
	h∇h =   [-x,-y,0]
	s² = x^2 + y^2
	J = jacobi(2s²-1,big(n)-1,big(3)//2, big(abs(m))//1 )

	ψp =   J * ((m < 0) ? CSH.sinsinpoly(-m, x, y) : CSH.cossinpoly(m, x, y))
	ψp = convert_polynom(T,ψp/sqrt(prefac(big(n),big(abs(m)))))

	return h²*∇(ψp)×ez + 3*ψp*h∇h×ez - z*∇(ψp)×h∇h
end

function NM(b::QGRIMBasis; m=-(N-1):(N-1))
	N = b.N
	ms = [m for m in m for n=1:(N-abs(m)-1)÷2]
	ns = [n for m in m for n=1:(N-abs(m)-1)÷2]
	return ns, ms
end


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



struct InsMFONCBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
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

const 𝐫 = [x, y, z]

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



## Basis following Ivers (2017) 

L(V::Ellipsoid) = [1/V.a 0 0
			0 1/V.b 0
			0 0 1/V.c]

L(V::Sphere) = I

radius(V::Sphere) = V.r₀
radius(V::Ellipsoid) = 1 #todo for the moment

const r² = 𝐫⋅𝐫

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

function t_mf_lapl(T,l,n) 
	# acobi(2r2-1, n, zero(T), T(l + 1//2))
    fac = 1/sqrt(l*(1 + l)*(1/T(-1 + 2*l + 4*n) + 1/T(3 + 2*l + 4*n)))
    return fac * (jacobi(2r2-1,n,T(0),T(l+1//2)) - jacobi(2r2-1,n-1,T(0),T(l+1//2))) 
end


function s_mf_lapl(T,l, n)
    # fac = (n==1) ? 1/sqrt(l*(1 + l)*(5 + 2*l)*(6 + l*(11 + 6*l))) : 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
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

#inviscid 3D basis, sphere.
function t_in(T,l,n) 
    fac = sqrt(T(3+2l+4n))/sqrt(T(l*(l+1)))
    return jacobi(2r2-1,n,T(0),T(l+1//2))*fac
end

function s_in(T,l,n) 
    fac = sqrt(T(5+2l+4n))/sqrt(T(4l*(l+1)*(n+1)^2))
    return (1-r2)*jacobi(2r2-1,n,T(1),T(l+1//2))*fac
end

utor_in(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (t_in(T,l,n)*rlm(l,m,x,y,z; norm=Full{T}, kwargs...)*𝐫)
upol_in(::Type{T}, l::Integer, m::Integer, n::Integer; kwargs...) where T = ∇ × (∇ × (s_in(T,l,n)*rlm(l,m,x,y,z; norm=Full{T}, kwargs...)*𝐫))


struct InviscidBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end


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
	lmnsp = [(l,m,n)  for l in 1:(N-1) for m in m for n in 0:((N-l+1)÷2-1) if abs(m)<=l]
	lmnst = [(l,m,n) for l in 1:N for m in m for n in 0:((N-l)÷2) if abs(m)<=l]
	return lmnsp, lmnst
end


## Generate constructors for each defined basis

for Basis in (:LebovitzBasis, :QGBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}; m = -N:N, kwargs...) where {T} =
                $Basis(N, V, basisvectors($Basis, N, V; kwargs...),false)
        ),
    )
end

for Basis in (:QGRIMBasis, :InsulatingMFBasis, :InsulatingMFCBasis, :GBasis, :PVMFBasis, :IversBasis, :IversBasisC, :InsMFLaplBasis, :InviscidBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} =
                $Basis(N, m, V, basisvectors($Basis, N, V; m, kwargs...),false)
        ),
    )
end

for Basis in (:QGIMBasis,:InsMFONBasis, :InsMFONCBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} =
                $Basis(N, m, V, basisvectors($Basis, N, V; m, kwargs...),true)
        ),
    )
end


isorthonormal(B::T) where T <: VectorBasis = B.orthonorm

## Geostrophic basis

# struct GBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
#     N::Int
#     V::Vol
#     el::Vector{vptype{T}}
#     orthonorm::Bool
# end


# function geo_veln(n::Integer,V::Volume{T}) where T

#     if typeof(V) <: Sphere
#         a, b, c = one(T), one(T), one(T)
#     elseif typeof(V) <: Ellipsoid
#         a, b, c = V.a, V.b, V.c
#     else
#         error("Geostrophic velocity not implemented for this Volume")
#     end

#     h2 = c^2*(1-x^2/a^2-y^2/b^2)
#     ez = [0,0,1]
#     hgradh = [-c^2*x/a^2,-c^2*y/b^2,0]
#     return (3+2n)//3*h2^n*z^0 * hgradh×ez
# end

# function basisvectors(::Type{GBasis}, N::Int, V::Volume{T}) where T
#     return [geo_veln(n,V) for n in 0:N÷2]
# end

