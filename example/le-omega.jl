using Distributed

addprocs()
@everywhere using BSON, Mire, TypedPolynomials, MultivariatePolynomials, LinearAlgebra, Statistics

datapath = "/cluster/home/gerickf/data/paperdata"

BSON.@load joinpath(datapath,"leomega_mire.bson") Les as us vs N

const n = length(as[1])
const a = 1.25
const b = 0.8
const c = 1.0

@everywhere begin
    MIREPATH=joinpath(dirname(pathof(Mire)),"..")
    include(joinpath(MIREPATH,"src/misc/analysis.jl"))
    include(joinpath(MIREPATH,"src/misc/plotting.jl"))


    function rmsp(us,N,vs,a,b,c,n,ngrid)
        [rmsle(N, vs, us[1:end÷2,i]/mean(us[:,i]), a, b, c,ngrid) for i=1:n]
    end
end

rms_all = pmap(u->rmsp(u,N,vs,a,b,c,n,30),us)
BSON.@save joinpath(datapath,"rms_all_mire.bson") rms_all
