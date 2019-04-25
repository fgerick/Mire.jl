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


function mat_force_galerkin!(A::AbstractArray{T,2},cmat,vs,N::Integer, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where T <: Real

    n_A = n_u(N)
    @assert size(A,1)==n_A
    @assert size(A,2)==n_A
    @assert length(vs)==n_A


    for j=1:n_A
        f = forcefun(vs[j],a,b,c,args...) #calculate f(uⱼ)
        for i=1:n_A
            # A[i,j] = inner_product(vs[i],f,a,b,c) # calculates ∫ <uᵢ,f(uⱼ)> dV
            A[i,j] = inner_product(cmat,vs[i],f; kwargs...)
        end
    end
end

function mat_force_galerkin_shift!(A::AbstractArray{T,2},cmat,vs,vs2,N::Int,Nshift::Int, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real

    n_A = n_u(N+Nshift)
    n_A2 = n_u(N)
    # @assert size(A,1)==n_A
    # @assert size(A,2)==n_A
    # @assert length(vs)==n_A


    for j=1:n_A2
        f = forcefun(vs2[j],a,b,c,args...) #calculate f(uⱼ)
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

function mat_force(N::Integer,cmat,vs, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where T <: Real
    n_combos = n_u(N)
    @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos)
    mat_force_galerkin!(A,cmat,vs,N ,forcefun,a,b,c,args...;kwargs...)
    return A
end
function mat_force_shift(N::Integer,Nshift::Integer,cmat,vs,vs2, forcefun::Function,a::T,b::T,c::T, args...) where T <: Real
    n_combos = n_u(N+Nshift)
    n_combos2 = n_u(N)

    # @assert n_combos == length(vs)
    A = spzeros(T,n_combos,n_combos2)
    mat_force_galerkin_shift!(A,cmat,vs,vs2,N,Nshift,forcefun,a,b,c,args...)
    return A
end

# """
#     assemblemhd(N,a,b,c,Ω,b0)
#
# Assembles MHD eigen system, such that
# λAx=Bx
#
# This is the dissipationless model, with
#
# ∂ₜu = -2Ω×u + (∇×b0)×b + (∇×b)×b0
# ∂ₜb = ∇×(u×b0)
#
# with Ω = 1/Le * eΩ.
# """
# function assemblemhd(N::Int,a::T,b::T,c::T,Ω,b0) where T <: Real
#     # T = typeof(a)
#     n_mat = n_u(N)
#     vs = vel(N,a,b,c)
#
#     A = spzeros(T,2n_mat,2n_mat)
#     B = spzeros(T,2n_mat,2n_mat)
#
#     A[1:n_mat,1:n_mat] .= mat_force(N,vs,inertial,a,b,c)'
#     A[n_mat+1:end,n_mat+1:end] .= mat_force(N,vs,inertialmag,a,b,c)'
#
#     B[1:n_mat,1:n_mat] .= mat_force(N,vs,coriolis,a,b,c,Ω)'
#     B[1:n_mat,n_mat+1:end] .= mat_force(N,vs,lorentz,a,b,c,b0)'
#
#     B[n_mat+1:end,1:n_mat] .= mat_force(N,vs,advection,a,b,c,b0)'
#
#     return A,B, vs
# end
function assemblemhd(N::Int,cmat,a::T,b::T,c::T,Ω,b0; kwargs...) where T<:Real
    # T = typeof(a)
    n_mat = n_u(N)
    vs = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat)
    B = spzeros(T,2n_mat,2n_mat)

    A[1:n_mat,1:n_mat] .= mat_force(N,cmat,vs,inertial,a,b,c; kwargs...)
    A[n_mat+1:end,n_mat+1:end] .= mat_force(N,cmat,vs,inertialmag,a,b,c; kwargs...)

    B[1:n_mat,1:n_mat] .= mat_force(N,cmat,vs,coriolis,a,b,c,Ω; kwargs...)
    B[1:n_mat,n_mat+1:end] .= mat_force(N,cmat,vs,lorentz,a,b,c,b0; kwargs...)

    B[n_mat+1:end,1:n_mat] .= mat_force(N,cmat,vs,advection,a,b,c,b0; kwargs...)

    return A,B, vs
end
function assemblemhd_shift(N::Int,Nshift::Int,cmat,a::T,b::T,c::T,Ω,b0; kwargs...) where T<:Real
    # T = typeof(a)
    n_mat = n_u(N+Nshift)
    n_mat2 = n_u(N)
    vs = vel(N+Nshift,a,b,c)
    vs2 = vel(N,a,b,c)

    A = spzeros(T,2n_mat,2n_mat2)
    B = spzeros(T,2n_mat,2n_mat2)

    A[1:n_mat,1:n_mat2] .= mat_force_shift(N,Nshift,cmat,vs,vs2,inertial,a,b,c; kwargs...)
    A[n_mat+1:end,n_mat2+1:end] .= mat_force_shift(N,Nshift,cmat,vs,vs2,inertialmag,a,b,c; kwargs...)

    B[1:n_mat,1:n_mat2] .= mat_force_shift(N,Nshift,cmat,vs,vs2,coriolis,a,b,c,Ω; kwargs...)
    B[1:n_mat,n_mat2+1:end] .= mat_force_shift(N,Nshift,cmat,vs,vs2,lorentz,a,b,c,b0; kwargs...)

    B[n_mat+1:end,1:n_mat2] .= mat_force_shift(N,Nshift,cmat,vs,vs2,advection,a,b,c,b0; kwargs...)

    return A,B, vs
end
#
# function assemblemhd_diffusion(N::Int,cmat,a::T,b::T,c::T,Ω,b0,Lu,Pm) where T<:Real
#     # T = typeof(a)
#     n_mat = n_u(N)
#     vs = vel(N,a,b,c)
#
#     A = spzeros(T,2n_mat,2n_mat)
#     B = spzeros(T,2n_mat,2n_mat)
#
#     A[1:n_mat,1:n_mat] .= mat_force(N,cmat,vs,inertial,a,b,c)
#     A[n_mat+1:end,n_mat+1:end] .= mat_force(N,cmat,vs,inertialmag,a,b,c)
#
#     B[1:n_mat,1:n_mat] .= mat_force(N,cmat,vs,coriolis,a,b,c,Ω) .+ mat_force(N,cmat,vs,viscous,a,b,c,Lu,Pm)
#     B[1:n_mat,n_mat+1:end] .= mat_force(N,cmat,vs,lorentz,a,b,c,b0)
#
#     B[n_mat+1:end,1:n_mat] .= mat_force(N,cmat,vs,advection,a,b,c,b0)
#     B[n_mat+1:end,n_mat+1:end] .= mat_force(N,cmat,vs,diffusion,a,b,c,b0,Lu)
#     return A,B, vs
# end

"""
    assemblehd(N,a,b,c,Ω,b0)

Assembles HD eigen system, such that
λAx=Bx

This is the inviscid model, with

∂ₜu = -2Ω×u.
"""
function assemblehd(N::Int,a::T,b::T,c::T,Ω) where T <: Real

    vs = vel(N,a,b,c)

    A = mat_force(N,vs,inertial,a,b,c)
    B = mat_force(N,vs,coriolis,a,b,c,Ω)
    return A,B, vs
end

function assemblehd(N::Int,cmat,a::T,b::T,c::T,Ω) where T <: Real

    vs = vel(N,a,b,c)

    A = mat_force(N,cmat,vs,inertial,a,b,c)
    B = mat_force(N,cmat,vs,coriolis,a,b,c,Ω)
    return A,B, vs
end
