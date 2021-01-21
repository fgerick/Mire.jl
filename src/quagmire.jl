
#Quagmire polynomial functions (2D reduced set of equations)

Π(n::Integer, m::Integer) = x^n*y^m
"""
    poly_inertial(n::Integer,m::Integer,a::T,b::T,c::T) where T

Polynomial for basis (n,m) of inertial force
"""
function poly_inertial(n::Integer,m::Integer,a::T,b::T,c::T) where T
    out = -((n+3)*(n+1)*c^2/a^2+(m+3)*(m+1)*c^2/b^2)*Π(n,m)+c^2*n*(n-1)*Π(n-2,m)
    out += c^2*m*(m-1)*Π(n,m-2)-n*(n-1)*c^2/b^2*Π(n-2,m+2)-m*(m-1)*c^2/a^2*Π(n+2,m-2)
    out += c^4/3*(n*(n-1)/b^4*Π(n-2,m+2) + m*(m-1)/a^4*Π(n+2,m-2)-(2n*m+n+m)/(a^2*b^2)*Π(n,m))
    return out
end

"""
    poly_coriolis(n::Integer,m::Integer,a::T,b::T,c::T,Ω::T) where T

Polynomial for basis (n,m) of Coriolis force.
"""
function poly_coriolis(n::Integer,m::Integer,a::T,b::T,c::T,Ω::T) where T
    E = -c^2/2*(x^2/a^2+y^2/b^2-1)
    ez = [0,0,1]
    return -2Ω*(∇(E) × ∇(Π(n,m)))⋅ez
end

#mhd:
"""
    poly_lorentz(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T

Polynomial for basis (n,m) of Lorentz force.
"""
function poly_lorentz(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T
    coeffs_A0, nm_pairs_A0 = coefficients(A0), exponents.(monomials(A0))
    D2A0 = sum([cf*poly_inertial(nmpair...,a,b,c) for (cf,nmpair) in zip(coeffs_A0,nm_pairs_A0)])
    D2A = poly_inertial(n,m,a,b,c)
    h2 = -c^2*(x^2/a^2+y^2/b^2-1)
    E = h2/2
    ez = [0,0,1]
    out =  h2.*(∇(D2A0)×∇(Π(n,m)))⋅ez .+ 3*Π(n,m) .*(∇(D2A0)×∇(E))⋅ez  .- D2A0.*(∇(E)×∇(Π(n,m)))⋅ez
    out += h2.*(∇(D2A) ×∇(A0))⋅ez     .+ 3*A0     .*(∇(D2A) ×∇(E))⋅ez  .- D2A .*(∇(E)×∇(A0))⋅ez
    return out
end

"""
    poly_maginertial(n::Integer,m::Integer,a::T,b::T,c::T) where T


"""
function poly_inertialmag(n::Integer,m::Integer,a::T,b::T,c::T) where T
    return Π(n,m)
end

"""
    poly_magadvection(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T


"""
function poly_advection(n::Integer,m::Integer,a::T,b::T,c::T,A0) where T
    h2 = -c^2*(x^2/a^2+y^2/b^2-1)
    E = h2/2
    ez = [0,0,1]
    out = h2.*(∇(Π(n,m)) × ∇(A0)) ⋅ ez .+ 3 .*A0 .*(∇(Π(n,m))×∇(E))⋅ez .+ 3 .*Π(n,m) .*(∇(E)×∇(A0))⋅ez
    return out
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
