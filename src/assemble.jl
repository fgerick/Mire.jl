### Functions to create matrices and assemble the full system.

## create matrices using galerkin:

"""
    mat_force_galerkin!(A,vs,N,forcefun,a,b,c, args...)

Fills Matrix `A` with Galerkin coefficients of force given by the function `forcefun(u,a,b,c,args...)`.
"""
function mat_force_galerkin!(A::AbstractArray{T,2},vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A
    @assert length(vs)==n_A

    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            A[i,j] = inner_product(vs[i],f,a,b,c) # calculates ∫ <uᵢ,f(uⱼ)> dV
        end
    end
end


function mat_force_galerkin!(A::AbstractArray{T,2},cmat,vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A
    @assert length(vs)==n_A


    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            # A[i,j] = inner_product(vs[i],f,a,b,c) # calculates ∫ <uᵢ,f(uⱼ)> dV
            A[i,j] = inner_product(cmat,vs[i],f)
        end
    end
end

"""
    mat_force(N,vs,forcefun,a,b,c, args...)

Allocates new matrix `A` and fills elements by calling
mat_force_galerkin!(A,vs,N,forcefun,a,b,c, args...).

Cached version:
mat_force(N,cmat,vs,forcefun,a,b,c, args...)

where `cmat[i,j,k]` contains the integrals of monomials xⁱyʲzᵏ.
"""
function mat_force(N::Integer,vs, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin!(A,vs,N ,forcefun,a,b,c,args...)
    return A
end

function mat_force(N::Integer,cmat,vs, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin!(A,cmat,vs,N ,forcefun,a,b,c,args...)
    return A
end


"""
    assemblemhd(N,a,b,c,Ω,b0)

Assembles MHD eigen system, such that
λAx=Bx

This is the dissipationless model, with

∂ₜu = -2Ω×u + (∇×b0)×b + (∇×b)×b0
∂ₜb = ∇×(u×b0)

with Ω = 1/Le * eΩ.
"""
function assemblemhd(N::Int,a::T,b::T,c::T,Ω,b0) where T <: Real
    # T = typeof(a)
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    A[1:n_mat,1:n_mat] .= mat_force(N,vs,inertial,a,b,c)
    A[n_mat+1:end,n_mat+1:end] .= mat_force(N,vs,inertial,a,b,c)

    B[1:n_mat,1:n_mat] .= mat_force(N,vs,coriolis,a,b,c,Ω)
    B[1:n_mat,n_mat+1:end] .= mat_force(N,vs,lorentz,a,b,c,b0)

    B[n_mat+1:end,1:n_mat] .= mat_force(N,vs,advection,a,b,c,b0)

    return A,B, vs
end
function assemblemhd(N::Int,cmat,a::T,b::T,c::T,Ω,b0) where T<:Real
    # T = typeof(a)
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    A[1:n_mat,1:n_mat] .= mat_force(N,cmat,vs,inertial,a,b,c)
    A[n_mat+1:end,n_mat+1:end] .= A[1:n_mat,1:n_mat] #mat_force_cached(N,cmat,vs,inertial,a,b,c)

    B[1:n_mat,1:n_mat] .= mat_force(N,cmat,vs,coriolis,a,b,c,Ω)
    B[1:n_mat,n_mat+1:end] .= mat_force(N,cmat,vs,lorentz,a,b,c,b0)

    B[n_mat+1:end,1:n_mat] .= mat_force(N,cmat,vs,advection,a,b,c,b0)

    return A,B, vs
end
