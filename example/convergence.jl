using Mire, BSON


datapath = "/cluster/home/gerickf/data/paperdata"

# datapath = "/Users/gerickf/data/paperdata/"


a,b,c = 1.,1.,1.
Le = 1e-6
Ω = 1/Le *Mire.ez
α = 0.01
v0_wu = [0,-xz/c^2,xy/b^2]
b0 = [-Mire.y/b^2,Mire.x/a^2,0] .+ α*v0_wu

N=3:2:35

for N in N
    cmat = Mire.cacheint(N,a,b,c)
    @time LHS,RHS, vs = assemblemhd(N,cmat,a,b,c,Ω,b0)
    BSON.@save joinpath(datapath,"mire_high_order_matrices_N$(N).bson") LHS RHS vs a b c Le
end
