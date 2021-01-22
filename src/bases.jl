
abstract type VectorBasis{T,V} end
abstract type Volume{T} end

ptype{T} =
    Polynomial{T,Term{T,Monomial{(x, y, z),3}},Array{Term{T,Monomial{(x, y, z),3}},1}}
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

Computes all combos i,j,k satisfying i+j+k < N.
"""
function combos(N::Int)
    [[i, j, k] for i = 0:N for j = 0:N for k = 0 if (i + j + k < N)],
    [[i, j, k] for i = 0:N for j = 0:N for k = 1:N if (i + j + k < N)]
end


"""
    basisvectors(::Type{LebovitzBasis},N::Int,V::Volume{T}) where T

Compute all velocity basis vectors for a given maximal degree N (Lebovitz 1989, eq. 41-42)
"""
function basisvectors(::Type{LebovitzBasis}, N::Int, V::Volume{T}) where {T}
    gp, hp = combos(N)
    v_1 = [v1(h..., V) for h in vcat(gp, hp)]
    v_2 = [v2(h..., V) for h in vcat(gp, hp)]
    v_3 = [v3(g..., V) for g in gp]

    return vcat(v_1, v_2, v_3)
end


## QG velocity basis
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



# Generate constructors for each defined basis
for Basis in (:LebovitzBasis, :QGBasis)
    eval(
        :(
            $Basis(N::Int, V::Volume{T}) where {T} =
                $Basis(N, V, basisvectors($Basis, N, V))
        ),
    )
end
