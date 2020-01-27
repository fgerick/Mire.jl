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
gammanp1half(n::Integer) = factorial(2n)//(4^n*factorial(n))

int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))


"""
    int_polynomial_ellipsoid(p, cmat)

DOCSTRING
"""
function int_polynomial_ellipsoid(p,cmat)
    ip = zero(eltype(cmat))
    cs = coefficients(p)
    exps = exponents.(monomial.(terms(p)))
    @simd for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end


"""
    inner_product(u,v,a,b,c)

Defines inner product in an ellipsoidal volume \$\\int\\langle u,v\\rangle dV\$.
"""
inner_product(u,v,a::Real,b::Real,c::Real) = int_polynomial_ellipsoid(dot(u,v),a,b,c)

"""
    inner_product(cmat, u, v; thresh=eps())

DOCSTRING

#Arguments:
- `cmat`: DESCRIPTION
- `u`: DESCRIPTION
- `v`: DESCRIPTION
- `thresh`: DESCRIPTION
"""
function inner_product(cmat, u, v; thresh=eps())
    if thresh != eps()
        u = truncpolyvec(u,thresh)
        v = truncpolyvec(v,thresh)
    end

    duv = dot(u,v)
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



# ∫xⁱyʲzᵏ(n×r)dS:

function int_surface_ellipsoid_torque(coordinate::Integer,i::Integer,j::Integer,k::Integer,a::T,b::T,c::T,r::T=one(T)) where T
        if coordinate==1
             if iseven(i) && isodd(i+j) && isodd(k)
                return 8*a^(1+i)*b^j*(c^2 - b^2)*c^k*r^(3 + i + j + k)*gamma((1+i)/2)*gamma(1+j/2)*gamma(1+k/2)/(4*gamma((5+i+j+k)/2))
            else
                return zero(T)
            end
        elseif coordinate==2
            if isodd(i) && isodd(i+j) && isodd(k)
                return 8*a^i*b^(1+j)*(a^2 - c^2)*c^k*r^(3 + i + j + k)*gamma(1+i/2)*gamma((1+j)/2)*gamma(1+k/2)/(4*gamma((5+i+j+k)/2))
            else
                return zero(T)
            end
        elseif coordinate==3
            if isodd(i) && isodd(j) && iseven(k)
                return 8*a^i*b^j*(b^2 - a^2)*c^(1 + k)*r^(3 + i + j + k)*gamma(1+i/2)*gamma(1+j/2)*gamma((1+k)/2)/(4*gamma((5+i+j+k)/2))
            else
                return zero(T)
            end
    else
        error("coordinate must be 1,2 or 3")
    end
end

function cacheint_surface_torque(N::Integer,coordinate::Integer,a::T,b::T,c::T,r::T=one(T)) where T
    return [int_surface_ellipsoid_torque(coordinate,big(i),big(j),big(k),a,b,c,r) for i=0:N,j=0:N,k=0:N]
end

cacheint_surface_torque(N::Integer,a::T,b::T,c::T,r::T = one(T)) where T = [cacheint_surface_torque(N,i,a,b,c,r) for i=1:3]
