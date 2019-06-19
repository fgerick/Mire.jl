function calcau(Les,N,cmat,a,b,c,Ω,b0)
    as = []
    us = []
    for Le in Les
        Ω = 1/Le * ez
        LHS,RHS,vs = Mire.assemblemhd(N,cmat,a,b,c,Ω,b0);
        esol = eigen(Matrix(RHS),Matrix(LHS));

        push!(as,esol.values)
        push!(us,esol.vectors)
    end
    return as,us

end

"""
proj(v,u)

projects v onto u, i.e projᵤv
"""
function proj(v,u)
    duv=dot(u,v)
    duu=dot(u,u)
    dd=duv/duu
    uout = [ui*dd for ui in u]
    return uout
end

function proj_without_norm(v,u)
    duv=dot(u,v)
    dd=duv
    uout = [ui*dd for ui in u]
    return uout
end


function geoageo(v,a,b,c)
    v_geo = proj(v,[-Mire.y/b^2,Mire.x/a^2,0])
    v_ageo = v - v_geo
    return v_geo, v_ageo
end

function rms_eqplane(v1,a,b,ngrid)
    X, Y = __mgrid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uz =  real.([v1[3](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    outsideellipse=radius.>1.
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    √sum(ux.^2+uy.^2+uz.^2)
end
function truncpoly(p,thresh=eps())
    c = coefficients(p)
    m = monomials(p)
    t = abs.(c) .> thresh
    return sum(c[t].*m[t])
end

truncpolyvec(v,thresh=eps()) = [truncpoly(vi,thresh=thresh) for vi in v]

function rms_ellipsoid(v1,a,b,c,ngrid)
    v = truncpolyvec(v1)
    X,Y,Z = meshgrid(range(-a,stop=a,length=ngrid), range(-b,stop=b,length=ngrid), range(-c,stop=c,length=ngrid))
    radius = .√(X.^2/a^2+Y.^2/b^2+Z.^2/c^2)
    insideellipse = radius .< 1.
    X = X[insideellipse]
    Y = Y[insideellipse]
    Z = Z[insideellipse]
    ux =  real.([v[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>zt) for (xt,yt,zt) in zip(X,Y,Z)])
    uy =  real.([v[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>zt) for (xt,yt,zt) in zip(X,Y,Z)])
    uz =  real.([v[3](Mire.x=>xt,Mire.y=>yt,Mire.z=>zt) for (xt,yt,zt) in zip(X,Y,Z)])
    √sum(ux.^2+uy.^2+uz.^2)
end

function rmsle(N,vs,u,a,b,c,ngrid=40)
    v_ile=Mire.eigenvel(N,vs,real.(u),a,b,c,norm=false)
    v_geo,v_ageo = geoageo(v_ile,a,b,c)
    rms2=rms_eqplane(v_geo,a,b,ngrid)
    rms3=rms_eqplane(v_ageo,a,b,ngrid)
    return Float64(rms3/rms2)
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
                  vz::AbstractVector{T}) where T
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = fill(1, m)
    on = fill(1, n)
    oo = fill(1, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end
