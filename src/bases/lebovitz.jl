

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

LebovitzBasis(N::Int, V::Volume{T}; m = -N:N, kwargs...) where {T} = LebovitzBasis(N, V, basisvectors(LebovitzBasis, N, V; kwargs...),false)