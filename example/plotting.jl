function grid(x,y)
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]
    return X,Y
end

function plot_velocity_equator(a,b,v1; ngrid=50, kwargs...)

    X, Y = grid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
#     uz =  real.([v1[3](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);

    outsideellipse=(X.^2/a^2+Y.^2/b^2).>0.99
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Y), Float64.(ux),
            Float64.(uy) ; color=Float64.(u), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_y(a,c,v1; ngrid=50, kwargs...)

    X, Z = grid(range(-a,stop=a,length=ngrid),range(-c,stop=c,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>0,Mire.z=>zt) for (xt,zt) in zip(X,Z)])
    uz =  real.([v1[3](Mire.x=>xt,Mire.y=>0,Mire.z=>zt) for (xt,zt) in zip(X,Z)])

    u = .√(ux.^2+uz.^2)

    outsideellipse=(X.^2/a^2+Z.^2/c^2).>=.99
    ux[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Z), Float64.(ux),
            Float64.(uz) ; color=Float64.(u), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_x(b,c,v1; ngrid=50, kwargs...)

    Y, Z = grid(range(-b,stop=b,length=ngrid),range(-c,stop=c,length=ngrid))
    uy =  real.([v1[2](Mire.x=>0,Mire.y=>yt,Mire.z=>zt) for (yt,zt) in zip(Y,Z)])
    uz =  real.([v1[3](Mire.x=>0,Mire.y=>yt,Mire.z=>zt) for (yt,zt) in zip(Y,Z)])

    u = .√(uy.^2+uz.^2)

    outsideellipse=(Y.^2/b^2+Z.^2/c^2).>=.99
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(Y), Float64.(Z), Float64.(uy),
            Float64.(uz) ; color=Float64.(u), kwargs...)


    ellipsex=b*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end


function plot_velocity_equator_sxyz(a,b,v1; ngrid=50, kwargs...)

    X, Y = grid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0,s=>sqrt(xt^2/a^2+yt^2/b^2)) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0,s=>sqrt(xt^2/a^2+yt^2/b^2)) for (xt,yt) in zip(X,Y)])
#     uz =  real.([v1[3](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);

    outsideellipse=(X.^2/a^2+Y.^2/b^2).>=1.0
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Y), Float64.(ux),
            Float64.(uy) ; color=Float64.(u), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_y_sxyz(a,c,v1; ngrid=50, kwargs...)

    X, Z = grid(range(-a,stop=a,length=ngrid),range(-c,stop=c,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>0,Mire.z=>zt,s=>xt/a) for (xt,zt) in zip(X,Z)])
    uz =  real.([v1[3](Mire.x=>xt,Mire.y=>0,Mire.z=>zt,s=>xt/a) for (xt,zt) in zip(X,Z)])

    u = .√(ux.^2+uz.^2)

    outsideellipse=(X.^2/a^2+Z.^2/c^2).>=.99
    ux[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Z), Float64.(ux),
            Float64.(uz) ; color=Float64.(u), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_x_sxyz(b,c,v1; ngrid=50, kwargs...)

    Y, Z = grid(range(-b,stop=b,length=ngrid),range(-c,stop=c,length=ngrid))
    uy =  real.([v1[2](Mire.x=>0,Mire.y=>yt,Mire.z=>zt,s=>yt/b) for (yt,zt) in zip(Y,Z)])
    uz =  real.([v1[3](Mire.x=>0,Mire.y=>yt,Mire.z=>zt,s=>yt/b) for (yt,zt) in zip(Y,Z)])

    u = .√(uy.^2+uz.^2)

    outsideellipse=(Y.^2/b^2+Z.^2/c^2).>=.99
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(Y), Float64.(Z), Float64.(uy),
            Float64.(uz) ; color=Float64.(u), kwargs...)


    ellipsex=b*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end


function plot_velocity_equator_Hsxyz(a,b,v1; ngrid=50, kwargs...)

    X, Y = grid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))

    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0,s=>sqrt(xt^2/a^2+yt^2/b^2),H=>sqrt(complex(1-xt^2/a^2-yt^2/b^2))) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0,s=>sqrt(xt^2/a^2+yt^2/b^2),H=>sqrt(complex(1-xt^2/a^2-yt^2/b^2))) for (xt,yt) in zip(X,Y)])
#     uz =  real.([v1[3](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);

    outsideellipse=(X.^2/a^2+Y.^2/b^2).>=1.0
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Y), Float64.(ux),
            Float64.(uy) ; color=Float64.(u), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_y_Hsxyz(a,c,v1; ngrid=50, kwargs...)

    X, Z = grid(range(-a,stop=a,length=ngrid),range(-c,stop=c,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>0,Mire.z=>zt,s=>xt/a,H=>sqrt(1-xt^2/a^2)) for (xt,zt) in zip(X,Z)])
    uz =  real.([v1[3](Mire.x=>xt,Mire.y=>0,Mire.z=>zt,s=>xt/a,H=>sqrt(1-xt^2/a^2)) for (xt,zt) in zip(X,Z)])

    u = .√(ux.^2+uz.^2)

    outsideellipse=(X.^2/a^2+Z.^2/c^2).>=.99
    ux[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(X), Float64.(Z), Float64.(ux),
            Float64.(uz) ; color=Float64.(u), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_x_Hsxyz(b,c,v1; ngrid=50, kwargs...)

    Y, Z = grid(range(-b,stop=b,length=ngrid),range(-c,stop=c,length=ngrid))
    uy =  real.([v1[2](Mire.x=>0,Mire.y=>yt,Mire.z=>zt,s=>yt/b,H=>sqrt(1-yt^2/b^2)) for (yt,zt) in zip(Y,Z)])
    uz =  real.([v1[3](Mire.x=>0,Mire.y=>yt,Mire.z=>zt,s=>yt/b,H=>sqrt(1-yt^2/b^2)) for (yt,zt) in zip(Y,Z)])

    u = .√(uy.^2+uz.^2)

    outsideellipse=(Y.^2/b^2+Z.^2/c^2).>=.99
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    streamplot(Float64.(Y), Float64.(Z), Float64.(uy),
            Float64.(uz) ; color=Float64.(u), kwargs...)


    ellipsex=b*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end
