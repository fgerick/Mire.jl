### Functions to define integrals of monomials over ellipsoidal volume and surface


"""
Integral over the surface of an ellipsoid.
"""
function int_ellipsoid_surface(p::Monomial,a::Real,b::Real,c::Real)
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    if iseven(i) && iseven(j) && iseven(k)

        f1 = factorial(i)*factorial(j)*factorial(k)
        f2 = factorial(i÷2)*factorial(j÷2)*factorial(k÷2)
        f3 = factorial(2+i+j+k)
        f4= factorial(1+(i+j+k)÷2)
        ft = f1/f2/f3*f4
        return √big(π)*a^(i+1)*b^(j+1)*c^(k+1)*ft
    else
        zero(BigFloat)
    end
end

"""
    int_monomial_ellipsoid(p,a,b,c)

Integrate monomial `p=xⁱyʲzᵏ` over ellipsoid of semi-axes `a,b,c`.
"""
function int_monomial_ellipsoid(p::Monomial,a::Real,b::Real,c::Real)
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    # i = exponent(p,x)
    # j = exponent(p,y)
    # k = exponent(p,z)
    return int_monomial_ellipsoid(i,j,k,a,b,c)
end

"""
    int_monomial_ellipsoid(i,j,k,a,b,c)

Integrate monomial `xⁱyʲzᵏ` over ellipsoid of semi-axes `a,b,c`.
"""
function int_monomial_ellipsoid(i::BigInt,j::BigInt,k::BigInt,a::Real,b::Real,c::Real)
    if iseven(i) && iseven(j) && iseven(k)
        γ₁ = i÷2
        γ₂ = j÷2
        γ₃ = k÷2
        γ = γ₁ + γ₂ + γ₃

        # : as stated in Jeremie's thesis and Vidal & Cebron: (not working)
        # f1 = factorial(γ+1)
        # f2 = factorial(2γ)
        # f3 = factorial(2γ+3)
        # f4 = factorial(γ₁)*factorial(γ₂)*factorial(γ₃)

        f1 = factorial(big(γ₁+γ₂+γ₃+1))
        f2 = factorial(2γ₁)*factorial(2γ₂)*factorial(2γ₃)
        f3 = factorial(2γ₁+2γ₂+2γ₃+3)
        f4 = factorial(γ₁)*factorial(γ₂)*factorial(γ₃)

        fact = f1*f2/f3/f4
        return 8big(π)*a^(2*γ₁+1)*b^(2*γ₂+1)*c^(2γ₃+1)*fact
    else
        zero(BigFloat)
    end
end


# function int_monomial_ellipsoid(i,j,k,a,b,c)
#     return ((1 + (-1)^i) *(1 + (-1)^j) *(1 + (-1)^k) *a^(1+i) *b^(1+j) *c^(1+k) *gamma((1 + i)/2)* gamma((1 + j)/2) * gamma((1 + k)/2))/
#                                                                         (8*gamma(1/2* (5 + i + j + k)))
# end



int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))
function int_polynomial_ellipsoid(p,cmat)
        ip = zero(eltype(cmat))
        cs = coefficients(p)
        exps = exponents.(monomial.(terms(p)))
        @inbounds @simd for i=1:length(cs)
            ip+=cs[i]*cmat[(exps[i] .+ 1)...]
        end
        return ip
    end
int_polynomial_ellipsoid_surface(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_ellipsoid_surface.(monomial.(terms(p)),a,b,c))

"""
    inner_product(u,v,a,b,c)

Defines inner product in an ellipsoidal volume \$\\int\\langle u,v\\rangle dV\$.
"""
inner_product(u,v,a::Real,b::Real,c::Real) = int_polynomial_ellipsoid(dot(u,v),a,b,c)

function inner_product(cmat,u,v)
    duv = dot(u,v)
    ip = zero(eltype(cmat))
    cs = coefficients(duv)
    exps = exponents.(monomial.(terms(duv)))
    @inbounds @simd for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end

"""
Function to precalculate monomial integrations.
"""
function cacheint(n::Int,a::T,b::T,c::T) where T<:Real
    Nmax=4n
    cachedmat=zeros(T,Nmax+1,Nmax+1,Nmax+1)
    for i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid(big(i),big(j),big(k),a,b,c)
    end
    return cachedmat
end

### Truncated integral

function int_monomial_ellipsoid_truncated(i::BigInt,j::BigInt,k::BigInt,a::Real,b::Real,c::Real,d::Real)
    if iseven(i) && iseven(j) && iseven(k)
            return a*b*c*d^3*(a*d)^i*(b*d)^j*(c*d)^k*gamma((1 + i)/2)*gamma((1 + j)/2) *gamma((1 + k)/2)/(8gamma(1/2* (5 + i + j + k)))
        else
        zero(BigFloat)
    end
end

function cacheint_truncated(n::Int,a::T,b::T,c::T,d::T) where T<:Real
    Nmax=4n
    cachedmat=zeros(T,Nmax+1,Nmax+1,Nmax+1)
    for i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid_truncated(big(i),big(j),big(k),a,b,c,d)
    end
    return cachedmat
end
