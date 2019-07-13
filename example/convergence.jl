using Mire, BSON


datapath = "/cluster/home/gerickf/data/paperdata"

# datapath = "/Users/gerickf/data/paperdata/"


a,b,c = 1.25,0.8,1.
Le = 1e-7
Ω = 1/Le *Mire.ez
α = 0.01
b0 = [-Mire.y/b^2,Mire.x/a^2,0] .+ α*[x^2,0,0]

N=3:2:35

for N in N
    cmat = Mire.cacheint(N,a,b,c)
    @time LHS,RHS, vs = assemblemhd(N,cmat,a,b,c,Ω,b0)
    BSON.@save joinpath(datapath,"mire_high_order_matrices_N$(N).bson") LHS RHS vs a b c Le
end
