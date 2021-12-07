module Quagmire
#Quagmire polynomial functions (2D reduced set of equations)

using LinearAlgebra
using MultivariatePolynomials
using TypedPolynomials
using SpecialFunctions
using SparseArrays

using Mire: ∇,Δ, x,y,z
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

## QG tools
qg_combos(N::Integer) = [[i,j] for i=0:N for j=0:N if (i+j<N)]
# qg_combos(N::Integer) = (N==0) ? [[0,0]] : vcat(combos(N-1),[[i,j] for i=0:N for j=0:N if (N-1<i+j<=N)]) #sorted

"""
    qg_vel(n::Integer,m::Integer,a::T,b::T,c::T) where T
"""
function qg_vel(n::Integer,m::Integer,a::T,b::T,c::T) where T
    h2 = c^2*(1-x^2/a^2-y^2/b^2)
    ez = [0,0,1]
    hgradh = [-c^2*x/a^2,-c^2*y/b^2,0]
    return h2*∇(x^n*y^m)×ez+3*x^n*y^m*hgradh×ez-z*∇(x^n*y^m)×hgradh
end

"""
    qg_vel(N::Integer,a::T,b::T,c::T) where T
Generate QG basis vectors.
"""
function qg_vel(N::Int,a::T,b::T,c::T) where T
    cs = qg_combos(N)
    return [qg_vel(ci...,a,b,c) for ci in cs]
end


#Quagmire (2D reduced equations)

function assemblemhd_quag(N::Int, a::T, b::T, c::T, Ω::T, A0;
    cmat = cacheint2D(N,a,b,c)) where T
vs_qg = qg_vel(N, a, b, c)
n_mat_qg = length(vs_qg)

nmat = 2n_mat_qg
A = spzeros(T, nmat, nmat)
B = spzeros(T, nmat, nmat)
projectforce_2D!(view(A, 1:n_mat_qg, 1:n_mat_qg),           N, cmat, poly_inertial,a,b,c)
projectforce_2D!(view(A, n_mat_qg+1:nmat, n_mat_qg+1:nmat), N, cmat, poly_inertialmag,a,b,c)
projectforce_2D!(view(B, 1:n_mat_qg, 1:n_mat_qg),           N, cmat, poly_coriolis,a,b,c, Ω)
projectforce_2D!(view(B, 1:n_mat_qg, n_mat_qg+1:nmat),      N, cmat, poly_lorentz,a,b,c, A0)
projectforce_2D!(view(B, n_mat_qg+1:nmat, 1:n_mat_qg),      N, cmat, poly_advection,a,b,c, A0)

return A, B, vs_qg
end



function projectforce_2D!(A::AbstractArray{T,2}, N::Integer, cmat, polyfun::Function,a::T,b::T,c::T, args...) where T

    combos = qg_combos(N)

    n_A = length(combos)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A

    @inbounds for j=1:n_A
        p = polyfun(combos[j]...,a,b,c,args...)
        for i=1:n_A
            A[i,j] = inner_product_2D(cmat,p,Π(combos[i]...),a,b,c)
        end
    end
    return nothing
end


function projectforce_2D(N::Integer,cmat, polyfun::Function, a::T,b::T,c::T, args...) where T
    n_combos = n_unknown(N)
    A = spzeros(T,n_combos,n_combos)
    projectforce_2D!(A,N,cmat, polyfun,a,b,c,args...)
    return A
end






end #module