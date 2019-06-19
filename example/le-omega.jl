using Distributed

addprocs()
@everywhere using BSON, Mire, TypedPolynomials, MultivariatePolynomials, LinearAlgebra, Statistics

datapath = "/cluster/home/gerickf/data/paperdata"

BSON.@load joinpath(datapath,"leomega_mire.bson") Les as us vs N

const n = length(as[1])
const a = 1.25
const b = 0.8
const c = 1.0

cmat = Mire.cacheint(N,a,b,c)

@everywhere begin
    MIREPATH=joinpath(dirname(pathof(Mire)),"..")
    include(joinpath(MIREPATH,"src/misc/analysis.jl"))
    include(joinpath(MIREPATH,"src/misc/plotting.jl"))


    function rmsp(us,N,vs,a,b,c,n,ngrid)
        [rmsle(N, vs, us[1:end÷2,i]/mean(us[:,i]), a, b, c,ngrid) for i=1:n]
    end

    function rms_parperp(N,vs,cmat,a,b,c,u)
        v_ile=Mire.eigenvel(N,vs,real.(u),a,b,c,norm=false)
        v_perp = v_ile[1]^2+v_ile[2]^2
        v_par = v_ile[3]^2
        rms_par=Mire.int_polynomial_ellipsoid(v_perp,cmat)
        rms_perp=Mire.int_polynomial_ellipsoid(v_par,cmat)
        return real(rms_perp)/real(rms_par)
    end

    function rmspp(N,vs,cmat,a,b,c,us,n)
        [rms_parperp(N, vs,cmat, a, b, c, us[1:end÷2,i]) for i=1:n]
    end
end

# rms_all = pmap(u->rmsp(u,N,vs,a,b,c,n,30),us)
rms_all = pmap(u->rmspp(N,vs,cmat,a,b,c,u,n),us)
BSON.@save joinpath(datapath,"rms_all_mire_parperp.bson") rms_all
