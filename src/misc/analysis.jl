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
