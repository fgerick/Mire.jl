function trackev_b_Mire(N,a,bs,c,α,Le,s0,s1,σ0,LHS0,RHS0,cmat; verbose=false, kwargs...)
    k=0
    bsout=[bs[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = []
    LHS = copy(LHS0)
    RHS = copy(RHS0)

    db = bs[2]-bs[1]
    dbt = db
    bt = bs[1]
    iter = 0
    # dlet = dle
    b_t = bs[1]

    while bsout[end] + db >= bs[end]

        if (k == 0 )
            b = bs[1]
        else
            b = bsout[end]+db
        end

        if b<bt
            b = bt
        end

        println()
        @show b
        flush(stdout)
        flush(stderr)
        cmat = Mire.cacheint(N,a,b,c)
        LHS,RHS = Mire.assemblemhd(N,cmat,1/b,b,c,1/Le * ez, B0_malkus(a,b));
        RHS[abs.(RHS).<1e-11].=0.
        LHS[abs.(LHS).<1e-11].=0.
        λ,u = eigstarget(RHS, LHS, target; v0 = utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            bt = bs[1] + db
        else
            corrs = [cor(real.(u[:,i]),real.(utarget)) for i=1:nev]
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort[1:3]]
            @show cors
                max_corr,imax = findmin(1.0 .- corrs)
            if corrs[imax] < 0.99 #if no correlating eigenvector is found the parameter stepping is too high
                @warn "Correlation is only $(corrs[imax])!, lowering step"
                dbt /= 2
                b_t = bsout[end] + dbt
                flush(stdout)
                flush(stderr)
                continue
            else
                push!(bsout,b)
                db = 2*(bsout[end]-bsout[end-1])
                b_t = bsout[end]+dbt
                dbt
                dbt = db
            end
        end
        push!(λs,λ[imax])
        push!(us,u[:,imax])
        target = λ[imax]
        utarget= u[:,imax]
        if verbose
            ω=abs(target)
            @show ω
            flush(stdout)
            flush(stderr)
        end
        push!(targets,target)
        k+=1
    end
    return λs,us,bsout
end
function trackev_b_Mire(N3,a,bs,c,α,Le,s0,s1,σ0,LHS0,RHS0,P0; verbose=false, kwargs...)
    k=0
    lesout=[les[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = []
    LHS = deepcopy(LHS0)
    RHS = deepcopy(RHS0)

    dle = les[1] #diff(les) # les[2]-les[1]
    dlet = dle
    le_t = les[1]
    iter = 0
    ms = -(P.nϕ-1)÷2:(P.nϕ-1)÷2
    mpeakold = 0
    while lesout[end] + dle <= les[end]

        if (k == 0 )
            Le = les[1]
        else
            Le = lesout[end]+dle #dle[k+1]
        end

        if Le>le_t
            Le = le_t
        end

        P = QG.MHDProblem(nh,nϕ,a,b,c,α,Le,P0.mean_field,s0,s1; approx_order=6)
            RHS[1:P.nh*P.nϕ,1:P.nh*nϕ] = QG.CoriolisForce(P.a, P.b, P.c, P.Le, P.H, P.nh, P.nϕ)
        λ,u = eigstarget(RHS, LHS, target; v0 = utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            le_t = les[1] + dle #[k+1]#[k+2]
            mpeakt = ms[findmax(sum(abs.(reshape(u[1:P.nh*P.nϕ,1],P.nh,P.nϕ)).^2,dims=1))[2][2]]
        else
            corrs = [cor(abs.(u[:,i]),abs.(utarget)) for i=1:nev]
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort]
            mpeaks = [ms[findmax(sum(abs.(reshape(u[1:P.nh*P.nϕ,i],P.nh,P.nϕ)).^2,dims=1))[2][2]] for i=1:nev]
                max_corr,imax = findmin(1.0 .- abs.(corrs))
            mpeakt=mpeaks[imax]
            if abs.(corrs[imax]) < 0.999 #if no correlating eigenvector is found the parameter stepping is too high
                if verbose
                    @warn "Correlation is only $(abs(corrs[imax]))!, lowering step"
                        end
                dlet /= 2
                le_t = lesout[end] + dlet
                flush(stdout)
                flush(stderr)
                continue
            else
                push!(lesout,Le)
                dle = 2*(lesout[end]-lesout[end-1])
                le_t = lesout[end] + dle #dle[k+1]
                dlet = dle #[k+1]
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        target = λ[imax]
        utarget= u[:,imax]
        mpeakold=mpeakt
        if verbose
            ω=abs(target)
            @show ω
            @show mpeakt
            flush(stdout)
            flush(stderr)
        end
        push!(targets,target)
        k+=1
        @show Le
        flush(stdout)
        flush(stderr)
    end
    return λs,us,lesout
end
