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
        coeff = gammanp1half(I)*gammanp1half(J)*gammanp1half(K)//(gammanp1half(D))
        return convert(T,a^(1 + i)*b^(1 + j)*c^(1 + k) *coeff)
    else
        return zero(T)
    end
end

#Γ(1/2+n)/√π
gammanp1half(n::Integer) = factorial(2n)//(big(4)^n*factorial(n))

# function int_ellipsoid(i,j,k,a::T,b::T,c::T) where T
#     if iseven(i) && iseven(j) && iseven(k)


#         return a*b*c*a^i*b^j*c^k*exp(lgamma((1+i)/2)+lgamma((1+j)/2)+lgamma((1+k)/2)-lgamma((5+i+j+k)/2))
#     else
#         return zero(T)
#     end
# end

# function int_monomial_ellipsoid2(i,j,k,a::T,b::T,c::T) where T
#     if iseven(i) && iseven(j) && iseven(k)

#         coeff = 32one(T)
#         for l in (4+i+j+k)÷2:(4+i+j+k)
#             coeff/=l
#         end
#         for l in (i,j,k)
#             if l!=0
#                 for l_ in l÷2:l
#                     coeff*=l_
#                 end
#             end
#         end

#         return coeff*a^(1 + i)*b^(1 + j)*c^(1 + k)
#     else
#         return zero(T)
#     end
# end
# function int_monomial_ellipsoid(i::Integer, j::Integer, k::Integer, a::AbstractFloat, b::AbstractFloat, c::AbstractFloat)
#     T = promote_type(typeof(a),typeof(b),typeof(c))
#     if iseven(i) && iseven(j) && iseven(k)
#         coeff = gamma(T(1 + i)/2)*gamma(T(1 + j)/2)*gamma(T(1 + k)/2)/gamma(T(5 + i + j + k)/2)
#         return convert(T,a^(1 + i)*b^(1 + j)*c^(1 + k)*coeff)
#     else
#         return zero(T)
#     end
# end


int_polynomial_ellipsoid(p,a::Real,b::Real,c::Real) = sum(coefficients(p).*int_monomial_ellipsoid.(monomial.(terms(p)),a,b,c))


"""
    int_polynomial_ellipsoid(p, cmat)

DOCSTRING
"""
function int_polynomial_ellipsoid(p::Polynomial{S},cmat::Array{T,3}) where {T,S}
    cs = coefficients(p)
    ip = zero(T)*zero(S)
    exps = exponents.(monomial.(terms(p)))

    @assert maxdegree(p)<=size(cmat,1)
    @inbounds for i=1:length(cs)
        ip+=cs[i]*cmat[(exps[i] .+ 1)...]
    end
    return ip
end



function inner_product(u, v, cmat)
    out = zero(eltype(cmat))*zero(coefficienttype(first(u)))*zero(coefficienttype(first(v)))
    @inbounds for (ui,vi) = zip(u,v)
        for ti in terms(ui), tj in terms(vi)
            m = monomial(ti)*monomial(tj)
            coeff = conj(coefficient(ti))*coefficient(tj)
            i_,j_,k_ = exponents(m)
            if iseven(i_) && iseven(j_) && iseven(k_) #only the even monomials have a nonzero integral
                out += coeff*cmat[i_+1,j_+1,k_+1]
            end
        end
    end

    return out
end

"""
    inner_product(u,v,a,b,c)

Defines inner product in an ellipsoidal volume \$\\int\\langle u,v\\rangle dV\$.
"""
function inner_product(u,v,a,b,c)

    out = zero(promote_type(typeof(a),typeof(b),typeof(c)))*zero(coefficienttype(first(u)))*zero(coefficienttype(first(v)))
    @inbounds for (ui,vi) = zip(u,v)
        for ti in terms(ui), tj in terms(vi)
            m = monomial(ti)*monomial(tj)
            coeff = conj(coefficient(ti))*coefficient(tj)
            i,j,k = exponents(m)
            if iseven(i) && iseven(j) && iseven(k) #only the even monomials have a nonzero integral
                out += coeff*int_monomial_ellipsoid(big(i),big(j),big(k),a,b,c)
            end
        end
    end

    return out
end





"""
    cacheintellipsoid(n::Int, a::T, b::T, c::T) where T

Function to precalculate monomial integrations in an ellipsoid (a,b,c).
NOTE: omits the factor π for convenience (has to be reintroduced if needed)
"""
function cacheintellipsoid(n::Int, a::T, b::T, c::T) where T
    Nmax = 4n
    cachedmat = zeros(T,Nmax+1,Nmax+1,Nmax+1)
    Threads.@threads for i = 0:Nmax
        for j = 0:Nmax, k = 0:Nmax
            cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid(big(i), big(j), big(k), a, b, c)
        end
    end
    return cachedmat
end

cacheint(n::Int, V::Ellipsoid{T}) where T = cacheintellipsoid(n,V.a,V.b,V.c)
cacheint(n::Int, V::Sphere{T}) where T = cacheintellipsoid(n,one(T),one(T),one(T))
