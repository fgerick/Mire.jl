using SymEngine, Mire, SpecialFunctions
#
# # revise()
# import Mire.cacheint, Mire.int_monomial_ellipsoid
#
# function cacheint(n::Int,a::T,b::T,c::T) where T
#     Nmax=4n
#     cachedmat=zeros(T,Nmax+1,Nmax+1,Nmax+1)
#     for i=0:Nmax,j=0:Nmax,k=0:Nmax
#         cachedmat[i+1,j+1,k+1] = int_monomial_ellipsoid(big(i),big(j),big(k),a,b,c)
#     end
#     return cachedmat
# end
#
#
# function int_monomial_ellipsoid(i::BigInt,j::BigInt,k::BigInt,a::T,b::T,c::T) where T
#     if iseven(i) && iseven(j) && iseven(k)
#         a^(1+i)*b^(1+j)*c^(1+k) *gamma((1 + i)/2)*gamma((1 + j)/2)*gamma((1 + k)/2)/(8*gamma((5+i+j+k)/2))
#     else
#         zero(BigFloat)
#     end
# end

revise()

n=5
@vars a b c ω
# a,b,c,ω = 1.,1.,1.,1e7
cmat=cacheint(n,a,b,c; dtype=BigFloat)

b0 = [-y/a^2,x/b^2,0]
Ω = [0,0,ω]
# revise()

@time LHS,RHS,vs = Mire.assemblemhd(n,cmat,a,b,c,Ω,b0)

LHS.colptr


sub(p)::Float64=subs(p,a=>1.25,b=>.8,c=>1.,ω=>1e6)
LHSn=broadcast(p->sub(p),LHS)
RHSn=broadcast(p->sub(p),RHS)

# LHSn=broadcast(P->Float64(subs(P,a=>1.,b=>1.,c=>1.,ω=>1e6)),LHS)
# RHSn=broadcast(P->Float64(subs(P,a=>1.,b=>1.,c=>1.,ω=>1e6)),RHS)


esol=eigen(inv(Matrix(LHSn))*RHSn)
