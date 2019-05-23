function track_lehnert(a,b,c,les,σ0,LHS0,RHS0,cmat; verbose=false, kwargs...)
    k=0
    lesout=[les[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = []
    LHS = deepcopy(LHS0)
    RHS = deepcopy(RHS0)

    dle = les[1]
    dlet = dle
    le_t = les[1]
    iter = 0
    mpeakold = 0
    while lesout[end] + dle <= les[end]

        if (k == 0 )
            Le = les[1]
        else
            Le = lesout[end]+dle
        end

        if Le>le_t
            Le = le_t
        end
        Ωvec = 1/Le*Mire.ez
        #only coriolis force depends on Le:
            RHS[1:end÷2,1:end÷2] = Mire.mat_force(N,cmat,vs,coriolis,a,b,c,Ωvec)
        λ,u = eigstarget(RHS, LHS, target; v0 = utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            le_t = les[1] + dle
        else
            corrs = [abs(cor(u[:,i],utarget)) for i =1:nev]
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort]
            max_corr,imax = findmin(1.0 .- abs.(corrs))

            if abs(corrs[imax]) < 0.98 #if no correlating eigenvector is found the parameter stepping is too high
                if verbose
                    @warn "Correlation is only $(abs(corrs[imax]))!, lowering step"
                flush(stdout)
                flush(stderr)
                end
                dlet /= 2
                le_t = lesout[end] + dlet
                continue
            else
                push!(lesout,Le)
                dle = 2*(lesout[end]-lesout[end-1])
                le_t = lesout[end] + dle
                dlet = dle
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
        if verbose
            @show Le
            flush(stdout)
            flush(stderr)
        end
    end
    return λs,us,lesout
end
