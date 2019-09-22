using SymEngine, LinearAlgebra, Mire, SpecialFunctions

n=5
@vars a b c ω
cmat=cacheint(n,a,b,c; dtype=BigFloat)

b0 = [-y/a^2,x/b^2,0]
Ω = [0,0,ω]

LHS,RHS,vs = Mire.assemblemhd(n,a,b,c,Ω,b0)


sub(p)::Float64=subs(p,a=>1.25,b=>.8,c=>1.,ω=>1e6)

LHSn=broadcast(p->sub(p),LHS)
RHSn=broadcast(p->sub(p),RHS)

esol=eigen(inv(Matrix(LHSn))*RHSn)
