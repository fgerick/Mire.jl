using Mire


function combos_qg_perp(N::Int)
    gpairs = Vector{Int}[]
    hpairs = Vector{Int}[]
    @inbounds for k=0, i=0:N-1, j=0:N-1
        if (i + j + k <= N-1)
            if k==0
                # if i==j==0
                #     push!(hpairs,[i,j,k])
                # end
                push!(gpairs,[i,j,k])
            else
                push!(hpairs,[i,j,k])
            end
        end
    end
    return gpairs, hpairs
end

function combos_qg_par(N::Int)
    gpairs = Vector{Int}[]
    hpairs = Vector{Int}[]
    @inbounds for k=[0,1], i=0:N-1, j=0:N-1
        if (i + j + k <= N-1)
            if k==0
                # if i==j==0
                #     push!(hpairs,[i,j,k])
                # end
                push!(gpairs,[i,j,k])
            else
                push!(hpairs,[i,j,k])
            end
        end
    end
    return gpairs, hpairs
end

F2(a,b,c) = -(Mire.x^2/a^2+Mire.y^2/b^2-1)
v1(n::Int,m::Int,l::Int,a,b,c) = Mire.∇(Mire.Π(n,m,l)*F2(a,b,c))×Mire.ex
v2(n::Int,m::Int,l::Int,a,b,c) = Mire.∇(Mire.Π(n,m,l)*F2(a,b,c))×Mire.ey
v3(n::Int,m::Int,l::Int,a,b,c) = Mire.∇(Mire.Π(n,m,l)*F2(a,b,c))×Mire.ez


function vel_qg(N::Int,a,b,c)
    # gpperp,hpperp = combos_qg_perp(N)

    gppar,hppar = combos_qg_par(N)
    v_1 = [v1(h...,a,b,c) for h in vcat(gppar,hppar)]
    v_2 = [v2(h...,a,b,c) for h in vcat(gppar,hppar)]
    v_3 = [v3(g...,a,b,c) for g in gppar]

    return vcat(v_1,v_2,v_3)
end

N=7
a,b,c=1.,1.,1.

vs = vel_qg(N,a,b,c)

length(vs)

function mat_force_galerkin_qg!(A::AbstractArray{T,2},cmat::Array{T,3},vs,
            N::Integer, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T <: Real}

    # n_A = n_u(N)
    # @assert size(A,1)==n_A
    # @assert size(A,2)==n_A
    # @assert length(vs)==n_A
    n_A = size(A,1)

    @inbounds for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            # A[i,j] = inner_product(vs[i],f,a,b,c) # calculates ∫ <uᵢ,f(uⱼ)> dV
            A[i,j] = inner_product(cmat,vs[i],f; kwargs...)
        end
    end
end



function mat_force_qg(N::Integer,cmat::Array{T,3},vs, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T <: Real}
    # n_combos = n_u(N)
    # @assert n_combos == length(vs)
    n_combos = length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin_qg!(A,cmat,vs,N ,forcefun,a,b,c,args...;kwargs...)
    return A
end

cmat = Mire.cacheint(N,a,b,c)
A=mat_force_qg(N,cmat,vs,Mire.inertial,a,b,c)


using Plots
theme(:juno)
spy(A)


function assemblemhd_qg(N::Int,cmat::Array{T,3},a::T,b::T,c::T,Ω,b0; kwargs...) where T<:Real
    # T = typeof(a)
    # n_mat = n_u(N)

    vs = vel_qg(N,a,b,c)
    n_mat = length(vs)
    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    # mat_force_galerkin!(A[1:n_mat,1:n_mat],cmat,vs,N ,inertial,a,b,c;kwargs...)
    # mat_force_galerkin!(A[n_mat+1:end,n_mat+1:end],cmat,vs,N ,inertialmag,a,b,c;kwargs...)
    #
    # mat_force_galerkin!(B[1:n_mat,1:n_mat],cmat,vs,N ,coriolis,a,b,c,Ω;kwargs...)
    # mat_force_galerkin!(B[1:n_mat,n_mat+1:end],cmat,vs,N ,lorentz,a,b,c,b0;kwargs...)
    # mat_force_galerkin!(B[n_mat+1:end,1:n_mat],cmat,vs,N ,advection,a,b,c,b0;kwargs...)

    A[1:n_mat,1:n_mat] .= mat_force_qg(N,cmat,vs,Mire.inertial,a,b,c; kwargs...)
    A[n_mat+1:end,n_mat+1:end] .= mat_force_qg(N,cmat,vs,Mire.inertialmag,a,b,c; kwargs...)

    B[1:n_mat,1:n_mat] .= mat_force_qg(N,cmat,vs,Mire.coriolis,a,b,c,Ω; kwargs...)
    B[1:n_mat,n_mat+1:end] .= mat_force_qg(N,cmat,vs,Mire.lorentz,a,b,c,b0; kwargs...)

    B[n_mat+1:end,1:n_mat] .= mat_force_qg(N,cmat,vs,Mire.advection,a,b,c,b0; kwargs...)

    return A,B, vs
end

N=3
Ω = 1/1e-7*[0,0,1]
b0 = [-Mire.y/b^2,Mire.x/a^2,0]
A,B,vs = assemblemhd_qg(N,cmat,a,b,c,Ω,b0);



esol = eigen(inv(Matrix(A))*B)

plot(sort(abs.(esol.values).+eps()),yscale=:log10,marker=1)

sort(abs.(esol.values)[abs.(esol.values).<1e7],rev=true)
