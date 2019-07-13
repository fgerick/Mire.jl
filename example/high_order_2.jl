using Distributed, DistributedArrays

addprocs(4)
@everywhere using BSON, Mire, SparseArrays, TypedPolynomials, MultivariatePolynomials, LinearAlgebra, Statistics

# datapath = "/cluster/home/gerickf/data/paperdata"

# N=parse(Int64,ARGS[1])
N=7
# const n = length(as[1])
const a = 1.25
const b = 0.8
const c = 1.0
const Le = 1e-7
const Ω = 1/Le *Mire.ez
const b0 = [-Mire.y/b^2,Mire.x/a^2,0]
# const N = 7
cmat = Mire.cacheint(N,a,b,c)

@everywhere begin
    # function galerkin_ij(i,j,vs,forcefun,a,b,c,cmat,args...;kwargs...)
    #     f = forcefun(vs[j],a,b,c,args...)
    #     inner_product(cmat,vs[i],f; kwargs...)
    # end
    # function mat_force_galerkin_parallel(cmat::Array{T,3},vs::Array{Array{P,1},1},
    #             N::Integer, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T <: Real, P <: Polynomial{T}}
    #
    #     n_A = n_u(N)
    #     # @assert size(A,1)==n_A
    #     # @assert size(A,2)==n_A
    #     # @assert length(vs)==n_A
    #     inds = [(i,j) for i=1:n_A,j=1:n_A]
    #     A = pmap(ind->galerkin_ij(ind...,vs,forcefun,a,b,c,cmat,args...;kwargs...),inds)
    #
    # end


    function assemblemhd_parallel(N::Int,cmat::Array{T,3},a::T,b::T,c::T,Ω,b0; kwargs...) where T<:Real
        # T = typeof(a)
        n_mat = n_u(N)
        vs = vel(N,a,b,c)

        A = spzeros(T,2n_mat,2n_mat)
        B = spzeros(T,2n_mat,2n_mat)


        A[1:n_mat,1:n_mat] .= mat_force_galerkin_parallel(cmat,vs,N,Mire.inertial,a,b,c; kwargs...)
        A[n_mat+1:end,n_mat+1:end] .= mat_force_galerkin_parallel(cmat,vs,N,Mire.inertialmag,a,b,c; kwargs...)

        B[1:n_mat,1:n_mat] .= mat_force_galerkin_parallel(cmat,vs,N,Mire.coriolis,a,b,c,Ω; kwargs...)
        B[1:n_mat,n_mat+1:end] .= mat_force_galerkin_parallel(cmat,vs,N,Mire.lorentz,a,b,c,b0; kwargs...)

        B[n_mat+1:end,1:n_mat] .= mat_force_galerkin_parallel(cmat,vs,N,Mire.advection,a,b,c,b0; kwargs...)

        return A,B, vs
    end

end

@time LHS,RHS,vs = assemblemhd_parallel(N,cmat,a,b,c,Ω,b0)

@time LHS2,RHS2, vs = assemblemhd(N,cmat,a,b,c,Ω,b0)

# BSON.@save joinpath(datapath,"high_order_N$(N).bson") LHS,RHS
