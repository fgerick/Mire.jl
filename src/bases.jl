
abstract type VectorBasis{T,V} end
abstract type Volume{T} end

ptype{T} = Polynomial{T,Term{T,Monomial{(x, y, z),3}},Array{Term{T,Monomial{(x, y, z),3}},1}}
vptype{T} = Vector{ptype{T}}

struct Ellipsoid{T<:Number} <: Volume{T}
    a::T
    b::T
    c::T
end

struct Sphere{T<:Number} <: Volume{T} end

Sphere() = Sphere{Float64}()
Ellipsoid(a, b, c) = Ellipsoid(promote(a, b, c)...)
Ellipsoid(a::Int, b::Int, c::Int) = Ellipsoid(a // 1, b // 1, c // 1)

struct LebovitzBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
end

struct QGBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
end

const ConductingMFBasis = LebovitzBasis

struct InsulatingMFBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}
    N::Int
    V::Vol
    el::Vector{vptype{T}}
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


"""
    combos(N::Int)

Computes all combos i,j,k satisfying i+j+k < N. Returns two `Vector{Vector{Int}}` with
the first one all `[i,j,k]` with `k=0`, and the second one with `1<=k<N`.
"""
function combos(N::Int)
    [[i, j, k] for i = 0:N for j = 0:N for k = 0 if (i + j + k < N)],
    [[i, j, k] for i = 0:N for j = 0:N for k = 1:N if (i + j + k < N)]
end


"""
    basisvectors(::Type{LebovitzBasis}, N::Int, V::Volume{T}) where T

Compute all velocity basis vectors for a given maximal degree `N` (Lebovitz 1989, eq. 41-42)
"""
function basisvectors(::Type{LebovitzBasis}, N::Int, V::Volume{T}) where {T}
    gp, hp = combos(N)
    v_1 = [v1(h..., V) for h in vcat(gp, hp)]
    v_2 = [v2(h..., V) for h in vcat(gp, hp)]
    v_3 = [v3(g..., V) for g in gp]

    return vcat(v_1, v_2, v_3)
end


## QG velocity basis (Gerick et al., 2020)

qg_combos(N::Integer) = [[i, j] for i = 0:N for j = 0:N if (i + j < N)]
# qg_combos(N::Integer) = (N==0) ? [[0,0]] : vcat(combos(N-1),[[i,j] for i=0:N for j=0:N if (N-1<i+j<=N)]) #sorted

"""
    uqg(n::Integer, m::Integer, V::Volume{T}) where T

Generate basis vector \$\\mathbf{u}=\\nabla(h^3x^ny^m)\\times\\nabla(z/h)\$
"""
function uqg(n::Integer, m::Integer, V::Volume{T}) where {T}

    if typeof(V) <: Sphere
        a, b, c = one(T), one(T), one(T)
    elseif typeof(V) <: Ellipsoid
        a, b, c = V.a, V.b, V.c
    else
        error("QG velocity not implemented for this Volume")
    end

    h2 = c^2 * (1 - x^2 / a^2 - y^2 / b^2)
    ez = [0, 0, 1]
    hgradh = [-c^2 * x / a^2, -c^2 * y / b^2, 0]
    return h2 * ∇(x^n * y^m) × ez + 3 * x^n * y^m * hgradh × ez - z * ∇(x^n * y^m) × hgradh
end


function basisvectors(::Type{QGBasis}, N::Int, V::Volume{T}) where {T}
    cs = qg_combos(N)
    return [uqg(ci..., V) for ci in cs]
end


## Insulating 3-D magnetic field basis in the sphere (Gerick et al., 2021)

#coefficients of interior and exterior potential field correction
vi(l::Int, n::Int) = -(l + 1)//(4l*n + 4l + 2n + 2)
ve(l::Int, n::Int) = l//((2*l + 1)*(2l + 2n + 3))

function basisvectors(::Type{InsulatingMFBasis}, N::Int, V::Volume{T}; norm=Schmidt{T}()) where T
    if typeof(V) != Sphere{T}
        return throw(ArgumentError("Insulating magnetic field basis is only implemented in the sphere!"))
    end
	r2 = x^2 + y^2 + z^2

	#l,m,n for poloidal field vectors
	ls = [l  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	ms = [m  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	ns = [n  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]

	NPOL = length(ls)

	#l,m,n for toroidal field vectors
	lstor = [l for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	mstor = [m for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	nstor  = [n for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]

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
    N = P.N
	#l,m,n for poloidal field vectors
	ls = [l  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	ms = [m  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]
	ns = [n  for l in 1:N for m in -N:N for n in 0:(N-l)÷2 if abs(m)<=l]

	#l,m,n for toroidal field vectors
	lstor = [l for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	mstor = [m for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]
	nstor  = [n for l in 1:(N-1) for m in -(N-1):(N-1) for n in 0:((N+1-l)÷2-1) if abs(m)<=l]

    return ls,ms,ns,lstor,mstor,nstor
end


# Generate constructors for each defined basis
for Basis in (:LebovitzBasis, :QGBasis, :InsulatingMFBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}; kwargs...) where {T} =
                $Basis(N, V, basisvectors($Basis, N, V; kwargs...))
        ),
    )
end
