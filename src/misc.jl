
#truncate polynomial
function truncpoly(u;atol=√(eps()))
    uo=zero(u)
    ti = terms(u)
    for t in ti
        c=coefficient(t)
        if abs(c)>atol
            uo+=t
        end
    end
    return uo
end

function truncvec(u;atol=√(eps()))
    uo=zero(u)
    for i=1:3
        ti = terms(u[i])
        for t in ti
            c=coefficient(t)
            if abs(c)>atol
                uo[i]+=t
            end
        end
    end
    return uo
end




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


#integrations for geostrophic velocity

function integrate_sxyzpoly(l,i,j,k,a,b,c)
    if iseven(k) && iseven(j) && iseven(i)
        return ((1 + (-1)^k)*a^(1 + i)*b^(1 + j)*c^(1 + k)*exp((im*i*π)/2)*(1 + exp(im*j*π))*π*
            gamma((1 + j)/2)*gamma((1 + k)/2)*gamma((2 + i + j + l)/2))/(4gamma((1-i)/2)*gamma((2 + i + j)/2)*gamma((5 + i + j + k + l)/2))
    else
        return zero(typeof(a))
    end
end

function cacheint_sxyzpoly(n::Int,a::T,b::T,c::T) where T
    Nmax=4n
    cachedmat=zeros(Complex{T},Nmax+1,Nmax+1,Nmax+1,Nmax+1)
    for l=0:Nmax,i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[l+1,i+1,j+1,k+1] = integrate_sxyzpoly(big(l),big(i),big(j),big(k),a,b,c)
    end
    return cachedmat
end


function integrate_Hsxyzpoly(m,l,i,j,k,a,b,c)
    if iseven(k) && iseven(j) && iseven(i)
        return ((1 + (-1)^k)*a^(1 + i)*b^(1 + j)*c^(1 + k)*exp((im*i*π)/2)*(1 + exp(im*j*π))*π*
            gamma((1 + j)/2)*gamma((2 + i + j + l)/2))*gamma((3 + k + m)/2)/(2*(1+k)*gamma((1-i)/2)*gamma((2 + i + j)/2)*gamma((5 + i + j + k + l + m)/2))
    else
        return zero(typeof(a))
    end
end

function cacheint_Hsxyzpoly(n::Int,a::T,b::T,c::T) where T
    Nmax=3n
    cachedmat=zeros(Complex{T},Nmax+1,Nmax+1,Nmax+1,Nmax+1,Nmax+1)
    @inbounds for m=0:Nmax,l=0:Nmax,i=0:Nmax,j=0:Nmax,k=0:Nmax
        cachedmat[m+1,l+1,i+1,j+1,k+1] = integrate_Hsxyzpoly(big(m),big(l),big(i),big(j),big(k),a,b,c)
    end
    return cachedmat
end
