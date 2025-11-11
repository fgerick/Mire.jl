

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


QGBasis(N::Int, V::Volume{T}; m = -N:N, kwargs...) where {T} = QGBasis(N, V, basisvectors(QGBasis, N, V; kwargs...),false)

"""
    QGIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}

Basis of complex QG inertial modes (Maffei et al., 2017). This basis is orthonormal.
"""
struct QGIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}
    N::Int
	m
    V::Vol
    el::Vector{vptype{T}}
    orthonorm::Bool
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

QGIMBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = QGIMBasis(N, m, V, basisvectors(QGIMBasis, N, V; m, kwargs...),true)

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


QGRIMBasis(N::Int, V::Volume{T}; m=-N:N, kwargs...) where {T} = QGRIMBasis(N, m, V, basisvectors(QGRIMBasis, N, V; m, kwargs...),false)