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

int_polynomial_ellipsoid_surface(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_ellipsoid_surface.(monomial.(terms(p)),a,b,c))

function int_monomial_ellipsoid(p::Monomial,a::Real,b::Real,c::Real)
    i = big(exponent(p,x))
    j = big(exponent(p,y))
    k = big(exponent(p,z))
    # i = exponent(p,x)
    # j = exponent(p,y)
    # k = exponent(p,z)
    return int_monomial_ellipsoid(i,j,k,a,b,c)
end

# """
#     int_monomial_ellipsoid(i,j,k,a,b,c)
#
# Integrate monomial `xⁱyʲzᵏ` over ellipsoid of semi-axes `a,b,c`.
# """
# function int_monomial_ellipsoid(i::BigInt,j::BigInt,k::BigInt,a::Real,b::Real,c::Real)
#     if iseven(i) && iseven(j) && iseven(k)
#         γ₁ = i÷2
#         γ₂ = j÷2
#         γ₃ = k÷2
#         γ = γ₁ + γ₂ + γ₃
#
#         # : as stated in Jeremie's thesis and Vidal & Cebron: (not working)
#         # f1 = factorial(γ+1)
#         # f2 = factorial(2γ)
#         # f3 = factorial(2γ+3)
#         # f4 = factorial(γ₁)*factorial(γ₂)*factorial(γ₃)
#
#         f1 = factorial(big(γ₁+γ₂+γ₃+1))
#         f2 = factorial(2γ₁)*factorial(2γ₂)*factorial(2γ₃)
#         f3 = factorial(2γ₁+2γ₂+2γ₃+3)
#         f4 = factorial(γ₁)*factorial(γ₂)*factorial(γ₃)
#
#         fact = f1*f2/f3/f4
#         return 8big(π)*a^(2*γ₁+1)*b^(2*γ₂+1)*c^(2γ₃+1)*fact
#     else
#         zero(BigFloat)
#     end
# end

"""
    int_monomial_ellipsoid(i::BigInt, j::BigInt, k::BigInt, a::T, b::T, c::T; dtype::DataType=BigFloat) where T

DOCSTRING

#Arguments:
- `i`: DESCRIPTION
- `j`: DESCRIPTION
- `k`: DESCRIPTION
- `a`: DESCRIPTION
- `b`: DESCRIPTION
- `c`: DESCRIPTION
- `dtype`: DESCRIPTION
"""
function int_monomial_ellipsoid(i::BigInt,j::BigInt,k::BigInt,a::T,b::T,c::T; dtype::DataType=Float64) where T
    if iseven(i) && iseven(j) && iseven(k)
        coeff=gamma((1 + i)/2)*gamma((1 + j)/2)*gamma((1 + k)/2)/(8*gamma((5+i+j+k)/2))
        (dtype == BigFloat) ? a^(1+i)*b^(1+j)*c^(1+k)*coeff : a^(1+i)*b^(1+j)*c^(1+k)*convert(dtype,coeff)
    else
        zero(dtype)
    end
end



int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))


"""
    int_polynomial_ellipsoid(p, cmat)

DOCSTRING
"""
function int_polynomial_ellipsoid(p,cmat)
    ip = zero(eltype(cmat))
    cs = coefficients(p)
    exps = exponents.(monomial.(terms(p)))
    @inbounds @simd for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end


#alternative dot product
function dotp(u,v)
    return u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
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
function inner_product(cmat,u,v; thresh=eps())
    if thresh != eps()
        u=truncpolyvec(u,thresh)
        v=truncpolyvec(v,thresh)
    end

    duv = dotp(u,v)
    return int_polynomial_ellipsoid(duv,cmat)
end


"""
Function to precalculate monomial integrations.
"""
function cacheint(n::Int,a::T,b::T,c::T; kwargs...) where T
    Nmax=4n
    cachedmat=zeros(T,Nmax+1,Nmax+1,Nmax+1)
    for i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid(big(i),big(j),big(k),a,b,c; kwargs...)
    end
    return cachedmat
end


#truncate polynomial

function truncpoly(p,thresh=eps())
    c = coefficients(p)
    m = monomials(p)
    t = abs.(c) .> thresh
    return sum(c[t].*m[t])
end

truncpolyvec(v,thresh=eps()) = [truncpoly(vi,thresh) for vi in v]



# function to integrate ∫x^iy^jz^k/(x^2/a^4+y^2/b^4)^2 dV
function int_monomial_ellipsoid_projected_ephi(i::BigInt,j::BigInt,k::BigInt,a::Real,b::Real,c::Real)
    if iseven(i) && iseven(j) && iseven(k) && ((i+j+k)>1)
        # a^(1+i)*b^(1+j)*c^(1+k) *gamma((1 + i)/2)*gamma((1 + j)/2)*gamma((1 + k)/2)/(4*(i+j)*gamma((3+i+j+k)/2))
        a^(1+i)*b^(1+j)*c^(1+k) *gamma((1 + i)/2)*gamma((1 + j)/2)*gamma((1 + k)/2)/((i+j-2)*(i+j)*(i+j+k-1)*gamma((i+j+k-1)/2))
    else
        zero(BigFloat)
    end
end

function cacheint_projected_ephi(n::Int,a::T,b::T,c::T) where T
    Nmax=4n
    cachedmat=zeros(T,Nmax+1,Nmax+1,Nmax+1)
    for i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid_projected_ephi(big(i),big(j),big(k),a,b,c)
    end
    return cachedmat
end
