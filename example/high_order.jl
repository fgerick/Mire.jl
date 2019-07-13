using Distributed

addprocs(24)
@everywhere using BSON, Mire, TypedPolynomials, MultivariatePolynomials, LinearAlgebra, Statistics

datapath = "/cluster/home/gerickf/data/paperdata"


const n = length(as[1])
const a = 1.25
const b = 0.8
const c = 1.0
const Le = 1e-7
const N = 15
cmat = Mire.cacheint(N,a,b,c)

@everywhere begin
    function galerkin_ij(vi,vj,forcefun,a,b,c,cmat,args...;kwargs...)
        f = forcefun(vj,a,b,c,args...)
        inner_product(cmat,vi,f; kwargs...)
    end
    function mat_force_galerkin_parallel(cmat::Array{T,3},vs::Array{Array{P,1},1},
                N::Integer, forcefun::Function,a::T,b::T,c::T, args...; kwargs...) where {T <: Real, P <: Polynomial{T}}

        n_A = n_u(N)
        @assert size(A,1)==n_A
        @assert size(A,2)==n_A
        @assert length(vs)==n_A

        A = pmap((i,j)->galerkin_ij(vs[i],vs[j],forcefun,a,b,c,cmat,args...;kwargs...) for i=1:n_A for j=1:n_A)

    end


    function assemblemhd_parallel(N::Int,cmat::Array{T,3},a::T,b::T,c::T,Ω,b0; kwargs...) where T<:Real
        # T = typeof(a)
        n_mat = n_u(N)
        vs = vel(N,a,b,c)

        A = spzeros(T,2n_mat,2n_mat)
        B = spzeros(T,2n_mat,2n_mat)


        A[1:n_mat,1:n_mat] .= mat_force_galerkin_parallel(cmat,vs,N,inertial,a,b,c; kwargs...)
        A[n_mat+1:end,n_mat+1:end] .= mat_force_galerkin_parallel(cmat,vs,N,inertialmag,a,b,c; kwargs...)

        B[1:n_mat,1:n_mat] .= mat_force_galerkin_parallel(cmat,vs,N,coriolis,a,b,c,Ω; kwargs...)
        B[1:n_mat,n_mat+1:end] .= mat_force_galerkin_parallel(cmat,vs,N,lorentz,a,b,c,b0; kwargs...)

        B[n_mat+1:end,1:n_mat] .= mat_force_galerkin_parallel(cmat,vs,N,advection,a,b,c,b0; kwargs...)

        return A,B, vs
    end

end

LHS,RHS,vs = assemblemhd_parallel(N,cmat,a,b,c,Ω,b0)

BSON.@save joinpath(datapath,"high_order_N$(N).bson") LHS,RHS
