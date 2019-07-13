using Distributed

addprocs()
@everywhere using BSON, Mire, TypedPolynomials, MultivariatePolynomials, LinearAlgebra, Statistics

datapath = "/cluster/home/gerickf/data/paperdata"

BSON.@load joinpath(datapath,"leomega_mire_N3.bson") Les as us vs N

const n = length(as[1])
const a = 1.25
const b = 0.8
const c = 1.0

cmat = Mire.cacheint(N,a,b,c)

@everywhere begin
    MIREPATH=joinpath(dirname(pathof(Mire)),"..")
    include(joinpath(MIREPATH,"src/misc/analysis.jl"))
    include(joinpath(MIREPATH,"src/misc/plotting.jl"))

    Car2Pol(x,y,a,b) = [√(x^2/a^2+y^2/b^2), atan(y/b,x/a)]
    function J(x,y,a,b,c)
        s,ϕ = Car2Pol(x,y,a,b)
        n1 = √(a^2*cos(ϕ)^2+b^2*sin(ϕ)^2)
        n2 = √(a^2*sin(ϕ)^2+b^2*cos(ϕ)^2)

         return [1/a*cos(ϕ)*n1 1/b*sin(ϕ)*n1 0
                 -1/a*sin(ϕ)*n2 1/b*cos(ϕ)*n2 0
                 0 0 1]
    end



function rmsle_2(N,vs,u,a,b,c,ngrid=40)
    v1=Mire.eigenvel(N,vs,real.(u/mean(u)),a,b,c,norm=false)

    ϕ = range(0,stop=2pi,length=ngrid)
    s = range(eps(),stop=1,length=ngrid)

    Φ,S = __mgrid(ϕ,s)
    X,Y = a*S.*cos.(Φ), b*S.*sin.(Φ)



    ux =  real.([v1[1](Mire.x=>xt, Mire.y=>yt, Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt, Mire.y=>yt, Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uz =  real.([v1[3](Mire.x=>xt, Mire.y=>yt, Mire.z=>0) for (xt,yt) in zip(X,Y)])

    Js = J.(X,Y,a,b,c)
    u_polar = [J*[ux,uy,uz] for (J,ux,uy,uz) in zip(Js,ux,uy,uz)];
    us,uphi,uz = [getindex.(u_polar,i) for i in 1:3];
    fact = @. sqrt(b^2*cos(Φ).^2 + a^2*sin(Φ)^2);
    uphi./=fact
    meanuphi = mean(uphi,dims=2)
    uphi .-= meanuphi
    uphi.*=fact

    rms_ug = √(ngrid*sum(meanuphi.^2))
#     rms_ag = .√(sum((uphi.-uphimean).^2 .+ us.^2 .+ uz.^2))
    rms_ag = √(sum(us.^2 .+ uphi.^2 .+ uz.^2))
    return rms_ag/rms_ug
end

    function rmsp(us,N,vs,a,b,c,n,ngrid)
        [rmsle(N, vs, us[1:end÷2,i], a, b, c,ngrid) for i=1:n]
    end

    # function rms_parperp(N,vs,cmat,a,b,c,u)
    #     v_ile=Mire.eigenvel(N,vs,real.(u),a,b,c,norm=false)
    #     v_perp = v_ile[1]^2+v_ile[2]^2
    #     v_par = v_ile[3]^2
    #     rms_par=Mire.int_polynomial_ellipsoid(v_perp,cmat)
    #     rms_perp=Mire.int_polynomial_ellipsoid(v_par,cmat)
    #     return real(rms_perp)/real(rms_par)
    # end
    #
    # function rmspp(N,vs,cmat,a,b,c,us,n)
    #     [rms_parperp(N, vs,cmat, a, b, c, us[1:end÷2,i]/mean(us[:,i])) for i=1:n]
    # end
end

rms_all = pmap(u->rmsp(u,N,vs,a,b,c,n,30),us)
# rms_all = pmap(u->rmspp(N,vs,cmat,a,b,c,u,n),us)
BSON.@save joinpath(datapath,"rms_all_mire_N$(N).bson") rms_all
