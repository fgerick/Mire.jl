
ptype{T} = Polynomial{T,Term{T,Monomial{(x, y, z),3}},Array{Term{T,Monomial{(x, y, z),3}},1}}
vptype{T} = Vector{ ptype{T}}

"""
    VectorBasis{T,V}

Defines `abstract type` of a vector basis to for a `T,V` so that `V=Volume{T}`.
"""
abstract type VectorBasis{T,V} end



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


function convert_polynom(S::Type{T},p) where T
	c = coefficients(p)
	m = monomials(p)
	return sum(S.(c).*m)
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

