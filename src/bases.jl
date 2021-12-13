
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
	ω = -1/(n*(2n+2m+1)+m/2+m^2/6)
	fac = gamma(n+3//2)*gamma(n+m)/(gamma(n)*(2n + 1//2 + m)*gamma(n + 3//2 + m))
	return -4big(pi)/ω*fac
end

function basiselementc(n::Integer,m::Integer,V::Sphere{T}) where T
	h² = 1-x^2-y^2
	ez =   [0,0,1]
	h∇h =   [-x,-y,0]
	s² = x^2 + y^2
	J = jacobi(2s²-1,big(n)-1,big(3)//2, big(m)//1 )

	ψp =   J * (im*CSH.sinsinpoly(big(m),x,y) + CSH.cossinpoly(big(m),x,y))
	ψp = convert_polynom(complex(T),ψp/sqrt(prefac(big(n),big(m))))

	return h²*∇(ψp)×ez + 3*ψp*h∇h×ez - z*∇(ψp)×h∇h
end

function basisvectors(::Type{QGIMBasis}, N::Int, V::Volume{T}) where T
    return [basiselementc(n,m,V) for m=0:(N-1) for n=1:(N-abs(m)-1)÷2 ]
end


function NM(b::QGIMBasis)
	N = b.N
	ms = [m for m=0:(N-1) for n=1:(N-abs(m)-1)÷2 ]
	ns = [n for m=0:(N-1) for n=1:(N-abs(m)-1)÷2 ]
	return ns,ms
end



"""
    QGRIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}

Real basis similar to the QG inertial modes (Maffei et al., 2017). They are not the solutions to the QG inertial mode equation and they are not orthogonal.
"""
struct QGRIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end


function basisvectors(::Type{QGRIMBasis}, N::Int, V::Volume{T}) where T
    return [basiselement(n,m,V) for m=-(N-1):(N-1) for n=1:(N-abs(m)-1)÷2]
end

function basiselement(n::Integer,m::Integer,V::Sphere{T}) where T
	h² = 1-x^2-y^2
	ez =   [0,0,1]
	h∇h =   [-x,-y,0]
	s² = x^2 + y^2
	J = jacobi(2s²-1,big(n)-1,big(3)//2, big(m)//1 )

	ψp =   J * ((m < 0) ? CSH.sinsinpoly(-m, x, y) : CSH.cossinpoly(m, x, y))
	ψp = convert_polynom(T,ψp/sqrt(prefac(big(n),big(abs(m)))))

	return h²*∇(ψp)×ez + 3*ψp*h∇h×ez - z*∇(ψp)×h∇h
end

function NM(b::QGRIMBasis)
	N = b.N
	ms = [m for m=-(N-1):(N-1) for n=1:(N-abs(m)-1)÷2]
	ns = [n for m=-(N-1):(N-1) for n=1:(N-abs(m)-1)÷2]
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
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end

#coefficients of interior and exterior potential field correction
vi(l::Int, n::Int) = -(l + 1)//(4l*n + 4l + 2n + 2)
ve(l::Int, n::Int) = l//((2*l + 1)*(2l + 2n + 3))

function bpol_g21(l,m,n; kwargs...)
	r² = x^2 + y^2 + z^2
	Rₗᵐ = rlm(l,m,x,y,z; kwargs...)
	Pₗₘₙ = -r²^(n+1)*Rₗᵐ / (2(n+1)*(2l+2n+3))
	return ∇×(∇×(Pₗₘₙ*[x,y,z])) - vi(l,n)*∇(Rₗᵐ)
end

function basisvectors(::Type{InsulatingMFBasis}, N::Int, V::Sphere{T}; norm=Schmidt{T}) where T
    if typeof(V) != Sphere{T}
        return throw(ArgumentError("Insulating magnetic field basis is only implemented in the sphere!"))
    end

	N-=1 # B of degree N instead of curl(B) of degree N

	r2 = x^2 + y^2 + z^2

	#l,m,n for poloidal field vectors
	lstor = [l  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	mstor = [m  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	nstor = [n  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]


	#l,m,n for toroidal field vectors
	ls = [l for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ms = [m for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ns  = [n for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]

	NPOL = length(ls)
	NTOR = length(lstor)

	#dummy variables to get datatypes for preallocation
	Qlm0 = r2^ns[1]*rlm(ls[1],ms[1],x,y,z; norm)
	Slm0 = (1-r2)*r2^nstor[1]*rlm(lstor[1],mstor[1],x,y,z; norm)
	Rlm0 = rlm(ls[1],ms[1],x,y,z; norm)*x^0*y^0*z^0

	j_tor0 = curl(Qlm0*[x,y,z])

	# preallocate vectors
	j_tor,b_pol = Vector{typeof(j_tor0)}(undef,NPOL),Vector{typeof(j_tor0)}(undef,NPOL)
	j_pol,b_tor = Vector{typeof(j_tor0)}(undef,NTOR),Vector{typeof(j_tor0)}(undef,NTOR)

	Qlm = zeros(typeof(Qlm0),NPOL)
	Slm = zeros(typeof(Slm0),NTOR)
	Rlm = zeros(typeof(Rlm0),NPOL)
	Plm = zeros(typeof(Qlm0),NPOL)

	#calculate all poloidal field vectors
	Threads.@threads for i = 1:NPOL
		n,m,l = ns[i],ms[i],ls[i]
	    Qlm[i] = r2^n*rlm(l,m,x,y,z,norm=norm)
		Rlm[i] = rlm(l,m,x,y,z,norm=norm)*x^0*y^0*z^0
		j_tor[i] = curl(Qlm[i]*[x,y,z])
		Plm[i] =  -r2^(n+1)*Rlm[i]/(2(n+1)*(2l+2n+3))
		b_pol[i] = curl(curl(Plm[i]*[x,y,z]))
		α = -vi(ls[i] ,ns[i])
		Vi = α*Rlm[i]
		b_pol[i] .+= ∇(Vi)
	end

	#calculate all toroidal field vectors
	@inbounds @simd for i = 1:NTOR
		Slm[i] = (1-r2)*r2^nstor[i]*rlm(lstor[i],mstor[i],x,y,z,norm=norm)
		b_tor[i] = curl(Slm[i]*[x,y,z])
		j_pol[i] = curl(b_tor[i])
	end

    # return j_pol, j_tor, b_pol, b_tor, Plm, Qlm, Slm, Rlm, ls, ms, ns, lstor, mstor, nstor
    return vcat(b_pol, b_tor)
end


function LMN(P::InsulatingMFBasis{T,V}) where {T,V}
	r2 = x^2 + y^2 + z^2
    N = P.N-1
	#l,m,n for poloidal field vectors

	#l,m,n for poloidal field vectors
	lstor = [l  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	mstor = [m  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	nstor = [n  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]


	#l,m,n for toroidal field vectors
	ls = [l for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ms = [m for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ns  = [n for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]


    return ls,ms,ns,lstor,mstor,nstor
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


function basisvectors(::Type{InsulatingMFCBasis}, N::Int, V::Sphere{T}; norm=Schmidt{T}) where T
	
	r2 = x^2 + y^2 + z^2

	N-=1
	#l,m,n for poloidal field vectors
	lstor = [l  for l in 1:N for m in 0:N for n in 0:(N-l)÷2 if abs(m)<=l]
	mstor = [m  for l in 1:N for m in 0:N for n in 0:(N-l)÷2 if abs(m)<=l]
	nstor = [n  for l in 1:N for m in 0:N for n in 0:(N-l)÷2 if abs(m)<=l]


	#l,m,n for toroidal field vectors
	ls = [l for l in 1:(N-1) for m in 0:(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ms = [m for l in 1:(N-1) for m in 0:(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ns  = [n for l in 1:(N-1) for m in 0:(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]


	NPOL = length(ls)
	NTOR = length(lstor)

	#dummy variables to get datatypes for preallocation
	Qlm0 = r2^ns[1]*rlm(ls[1],ms[1],x,y,z; norm,real=false)
	Slm0 = (1-r2)*r2^nstor[1]*rlm(lstor[1],mstor[1],x,y,z; norm, real=false)
	Rlm0 = rlm(ls[1],ms[1],x,y,z; norm, real=false)*x^0*y^0*z^0

	j_tor0 = curl(Qlm0*[x,y,z])

	# preallocate vectors
	j_tor,b_pol = Vector{typeof(j_tor0)}(undef,NPOL),Vector{typeof(j_tor0)}(undef,NPOL)
	j_pol,b_tor = Vector{typeof(j_tor0)}(undef,NTOR),Vector{typeof(j_tor0)}(undef,NTOR)

	Qlm = zeros(typeof(Qlm0),NPOL)
	Slm = zeros(typeof(Slm0),NTOR)
	Rlm = zeros(typeof(Rlm0),NPOL)
	Plm = zeros(typeof(Qlm0),NPOL)

	#calculate all poloidal field vectors
	Threads.@threads for i = 1:NPOL
		n,m,l = ns[i],ms[i],ls[i]
	    Qlm[i] = r2^n*rlm(l,m,x,y,z; norm=norm, real=false)
		Rlm[i] = rlm(l,m,x,y,z; norm=norm, real=false)*x^0*y^0*z^0
		j_tor[i] = curl(Qlm[i]*[x,y,z])
		Plm[i] =  -r2^(n+1)*Rlm[i]/(2(n+1)*(2l+2n+3))
		b_pol[i] = curl(curl(Plm[i]*[x,y,z]))
		α = -vi(ls[i] ,ns[i])
		Vi = α*Rlm[i]
		b_pol[i] .+= ∇(Vi)
	end

	#calculate all toroidal field vectors
	@inbounds @simd for i = 1:NTOR
		Slm[i] = (1-r2)*r2^nstor[i]*rlm(lstor[i],mstor[i],x,y,z; norm=norm,real=false)
		b_tor[i] = curl(Slm[i]*[x,y,z])
		j_pol[i] = curl(b_tor[i])
	end

    # return j_pol, j_tor, b_pol, b_tor, Plm, Qlm, Slm, Rlm, ls, ms, ns, lstor, mstor, nstor
    return vcat(b_pol, b_tor)
end

function LMN(P::InsulatingMFCBasis{T,V}) where {T,V}
	r2 = x^2 + y^2 + z^2
    N = P.N-1
	#l,m,n for poloidal field vectors
	lstor = [l  for l in 1:N for m in 0:N for n in 0:(N-l)÷2 if abs(m)<=l]
	mstor = [m  for l in 1:N for m in 0:N for n in 0:(N-l)÷2 if abs(m)<=l]
	nstor = [n  for l in 1:N for m in 0:N for n in 0:(N-l)÷2 if abs(m)<=l]


	#l,m,n for toroidal field vectors
	ls = [l for l in 1:(N-1) for m in 0:(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ms = [m for l in 1:(N-1) for m in 0:(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	ns  = [n for l in 1:(N-1) for m in 0:(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]

    return ls,ms,ns,lstor,mstor,nstor
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

btor(::Type{T}, n::Integer, m::Integer, l::Integer; kwargs...) where T = ∇ × (h(T,l,n)*rlm(l,m,x,y,z; norm=Schmidt{T}, kwargs...)*𝐫)
bpol(::Type{T}, n::Integer, m::Integer, l::Integer; kwargs...) where T = ∇ × (∇ × (k(T,l,n)*rlm(l,m,x,y,z; norm=Schmidt{T}, kwargs...)*𝐫))


function basisvectors(::Type{InsMFONCBasis}, N::Int, V::Sphere{T}; kwargs...) where T
	r2 = x^2+y^2+z^2
	N-=1

	ls = [l  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ms = [m  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ns = [n  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]

	NPOL = length(ls)

	lstor = [l for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	mstor = [m for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	nstor  = [n for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]

	NTOR = length(lstor)

	BP = map((n,l,m)->bpol(T,n,m,l; real=false, kwargs...),ns,ls,ms)
	BT = map((n,l,m)->btor(T,n,m,l; real=false, kwargs...),nstor,lstor,mstor)
    return vcat(BP,BT)
	# return BP, BT, ls, ms, ns, lstor,mstor,nstor
end



function LMN(P::InsMFONCBasis{T,V}) where {T,V}
	r2 = x^2 + y^2 + z^2
    N = P.N-1

	ls = [l  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ms = [m  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ns = [n  for l in 1:N for m in 0:N for n in 1:(N-l+2)÷2 if abs(m)<=l]

	lstor = [l for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	mstor = [m for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	nstor  = [n for l in 1:(N-1) for m in 0:(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]

    return ls,ms,ns,lstor,mstor,nstor
end



struct InsMFONBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
end


function basisvectors(::Type{InsMFONBasis}, N::Int, V::Sphere{T}; kwargs...) where T
	N-=1
	r2 = x^2+y^2+z^2
	ls = [l  for l in 1:N for m in -N:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ms = [m  for l in 1:N for m in -N:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ns = [n  for l in 1:N for m in -N:N for n in 1:(N-l+2)÷2 if abs(m)<=l]

	NPOL = length(ls)

	lstor = [l for l in 1:(N-1) for m in -(N-1):(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	mstor = [m for l in 1:(N-1) for m in -(N-1):(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	nstor  = [n for l in 1:(N-1) for m in -(N-1):(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]

	NTOR = length(lstor)

	BP = map((n,l,m)->bpol(T,n,m,l; kwargs...),ns,ls,ms)
	BT = map((n,l,m)->btor(T,n,m,l; kwargs...),nstor,lstor,mstor)
    return vcat(BP,BT)
	# return BP, BT, ls, ms, ns, lstor,mstor,nstor
end


function LMN(P::InsMFONBasis{T,V}) where {T,V}
	r2 = x^2 + y^2 + z^2
    N = P.N-1

	ls = [l  for l in 1:N for m in -N:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ms = [m  for l in 1:N for m in -N:N for n in 1:(N-l+2)÷2 if abs(m)<=l]
	ns = [n  for l in 1:N for m in -N:N for n in 1:(N-l+2)÷2 if abs(m)<=l]

	lstor = [l for l in 1:(N-1) for m in -(N-1):(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	mstor = [m for l in 1:(N-1) for m in -(N-1):(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]
	nstor  = [n for l in 1:(N-1) for m in -(N-1):(N-1) for n in 1:((N+1-l)÷2) if abs(m)<=l]

    return ls,ms,ns,lstor,mstor,nstor
end


# Generate constructors for each defined basis
for Basis in (:LebovitzBasis, :QGBasis, :QGRIMBasis, :InsulatingMFBasis, :InsulatingMFCBasis, :GBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}; kwargs...) where {T} =
                $Basis(N, V, basisvectors($Basis, N, V; kwargs...),false)
        ),
    )
end

for Basis in (:QGIMBasis,:InsMFONBasis, :InsMFONCBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}; kwargs...) where {T} =
                $Basis(N, V, basisvectors($Basis, N, V; kwargs...),true)
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

