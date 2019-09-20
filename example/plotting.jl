function __mgrid(x::AbstractVector{T},y::AbstractVector{T}) where T<: Real
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]
    return X,Y
end

function plot_velocity_equator_uphi(a,b,v1; ngrid=50, kwargs...)

    X, Y = __mgrid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);
    uphi =( -b*sin.(phi).*ux .+ a*cos.(phi).*uy)

    outsideellipse=(X.^2/a^2+Y.^2/b^2).>=1.
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    uphi[outsideellipse].=0.
    uphio = copy(uphi)

    # uphi.-=maximum(uphi)

    streamplot(Float64.(X), Float64.(Y), Float64.(ux),
            Float64.(uy) ; color=Float64.(uphi), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=1)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_equator_uz(a,b,v1; ngrid=50, kwargs...)

    X, Y = __mgrid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](x=>xt,y=>yt,z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](x=>xt,y=>yt,z=>0) for (xt,yt) in zip(X,Y)])
    uz =  real.([v1[3](x=>xt,y=>yt,z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);
    uphi =( -b*sin.(phi).*ux .+ a*cos.(phi).*uy)

    outsideellipse=radius.>=1.
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    uphi[outsideellipse].=0.
    uphio = copy(uphi)

    uphi.-=maximum(uphi)

    streamplot(Float64.(X), Float64.(Y), Float64.(ux), Float64.(uy) ; color=Float64.(uz), kwargs...) # linewidth=3*u/maximum(u), kwargs...)

    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=1)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_equator(a,b,v1; ngrid=50, kwargs...)

    X, Y = __mgrid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](x=>xt,y=>yt,z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](x=>xt,y=>yt,z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);
    uphi =( -b*sin.(phi).*ux .+ a*cos.(phi).*uy)

    outsideellipse=(X.^2/a^2+Y.^2/b^2).>=1.
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Y), Float64.(ux), Float64.(uy) ; color=Float64.(u), kwargs...) # linewidth=3*u/maximum(u), kwargs...)

    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=1)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end
function plot_velocity_at_z(a,b,zv,v1; ngrid=50, kwargs...)

    X, Y = __mgrid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](x=>xt,y=>yt,z=>zv) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](x=>xt,y=>yt,z=>zv) for (xt,yt) in zip(X,Y)])
    uz =  real.([v1[3](x=>xt,y=>yt,z=>zv) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);
    uphi =( -b*sin.(phi).*ux .+ a*cos.(phi).*uy)

    outsideellipse=radius.>=√(1.0-zv^2)
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    uphi[outsideellipse].=0.
    uphio = copy(uphi)

    uphi.-=maximum(uphi)

    streamplot(Float64.(X), Float64.(Y), Float64.(ux), Float64.(uy) ; color=Float64.(uz), kwargs...) # linewidth=3*u/maximum(u), kwargs...)

    ellipsex=a*√(1.0-zv^2)*cos.(range(0,stop=2π,length=100))
    ellipsey=b*√(1.0-zv^2)*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=1)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end
