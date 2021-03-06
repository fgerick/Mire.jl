### Functions to define integrals of monomials over ellipsoidal volume and surface

function int_monomial_ellipsoid(p::Monomial,a::T,b::T,c::T) where T
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    return int_monomial_ellipsoid(i,j,k,a,b,c)
end

"""
    int_monomial_ellipsoid(i::BigInt, j::BigInt, k::BigInt, a::T, b::T, c::T; dtype::DataType=BigFloat) where T

Integrate a monomial over the ellipsoids volume \$\\int x^iy^jz^k dV\$.

!!! warning "Missing factor"
     The factor \$\\pi\$ is removed in the integration due to numerical reasons.
     Remember to reintroduce it when the actual value of the integration is needed!

#Arguments:
- `i`: Exponent of `x`
- `j`: Exponent of `y`
- `k`: Exponent of `z`
- `a`: Semi-axis in `x`
- `b`: Semi-axis in `y`
- `c`: Semi-axis in `z`
"""
function int_monomial_ellipsoid(i::Integer,j::Integer,k::Integer,a::T,b::T,c::T) where T
    if iseven(i) && iseven(j) && iseven(k)
        I = i÷2
        J = j÷2
        K = k÷2
        D = (4+i+j+k)÷2
        coeff = gammanp1half(I)*gammanp1half(J)*gammanp1half(K)//(8*gammanp1half(D))
        return convert(T,(1 + (-1)^i)*(1 + (-1)^j)*(a^(1 + i))*(b^(1 + j))*c*((-c)^k + c^k) *coeff)
    else
        return zero(T)
    end
end

#Γ(1/2+n)/√π
gammanp1half(n::Integer) = factorial(2n)//(big(4)^n*factorial(n))

int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))


"""
    int_polynomial_ellipsoid(p, cmat)

DOCSTRING
"""
function int_polynomial_ellipsoid(p::Polynomial{S},cmat::Array{T,3}) where {T,S}
    ip = zero(promote_type(S,T))
    cs = coefficients(p)
    exps = exponents.(monomial.(terms(p)))

    @assert maxdegree(p)<=size(cmat,1)
    @inbounds for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end


"""
    inner_product(u,v,a,b,c)

Defines inner product in an ellipsoidal volume \$\\int\\langle u,v\\rangle dV\$.
"""
inner_product(u,v,a::Real,b::Real,c::Real) = int_polynomial_ellipsoid(dot(u,v),a,b,c)


dotp(u,v) = u[1]*v[1]+u[2]*v[2]+u[3]*v[3]

"""
    inner_product(cmat, u, v; thresh=eps())

DOCSTRING

#Arguments:
- `cmat`: DESCRIPTION
- `u`: DESCRIPTION
- `v`: DESCRIPTION
- `thresh`: DESCRIPTION
"""
function inner_product(cmat, u, v)
    duv = dot(u,v)
    return int_polynomial_ellipsoid(duv,cmat)
end

function inner_product_real(cmat, u, v)
    duv = dotp(u,v)
    return int_polynomial_ellipsoid(duv,cmat)
end

"""
Function to precalculate monomial integrations.
"""
function cacheint(n::Int, a::T, b::T, c::T) where T
    Nmax = 4n
    cachedmat = zeros(T,Nmax+1,Nmax+1,Nmax+1)
    @inbounds for i = 0:Nmax,j = 0:Nmax, k = 0:Nmax
        cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid(big(i), big(j), big(k), a, b, c)
    end
    return cachedmat
end

#Quagmire integration


# integrates ∫∫ x^i y^j sabc h^3 dsdϕ
function int_monomial_ellipse(i::BigInt,j::BigInt,a::Real,b::Real,c::Real)
    return ((3*(1+(-1)^i)*(1+(-1)^j)*a^(1+i)*b^(1+j)*c^3*√big(π)*gamma((1+i)/2)*gamma((1+j)/2))/(16*gamma((7+i+j)/2)))
end

function 𝒟(Π,a::T,b::T,c::T) where T
    if Π == 0
        return 0
    else
        cs = coefficients(Π)
        exps = exponents.(monomial.(terms(Π)))
        return sum([co*poly_inertial(nmpair...,a,b,c) for (nmpair,co) in zip(exps,cs)])
    end
end

function inner_product_2D(cmat,u,v,a::T,b::T,c::T) where T
    duv = u*𝒟(v,a,b,c)
    ip = zero(eltype(cmat))
    cs = coefficients(duv)
    exps = exponents.(monomials(duv))
    @inbounds for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end

function cacheint2D(n::Int,a::T,b::T,c::T) where T
    Nmax = 4n
    cmat = zeros(T,Nmax+1,Nmax+1)
    for i=0:Nmax,j=0:Nmax
        cmat[i+1,j+1] = int_monomial_ellipse(big(i),big(j),a,b,c)
    end
    return cmat
end
