#Solving for Malkus modes in the ellipsoid using a hybrid model

using Mire, LinearAlgebra

N,a,b,c = 7,1.25,0.8,1.0
cmat = cacheint(N,a,b,c)
vbasis = qg_vel(N,a,b,c)
bbasis = vel(N,a,b,c)

Le = 1e-4
Ω = [0,0,1/Le]
b0 = [-y/b^2,x/a^2,0] #modified Malkus field.

LHS,RHS = assemblemhd(Ω,b0,vbasis,bbasis,cmat)

#convert sparse to dense matrices and use LAPACK to solve generalized eigenproblem
esol = eigen(Matrix(RHS),Matrix(LHS))

#eigenvalues are complex with vanishingly small real part. Take only positive
#imaginary values:
evals = imag.(esol.values)[imag.(esol.values).>0]
evals = sort(evals,rev=true)

println("The ten largest eigen frequencies are:\n")
println(evals[1:10])

println("\n The ten smallest eigen frequencies are:\n")
println(reverse(evals[end-10:end]))
