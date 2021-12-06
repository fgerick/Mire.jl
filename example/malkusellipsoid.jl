#Solving for Malkus modes in the ellipsoid using a hybrid model

using Mire

N,a,b,c = 7,1.25,0.8,1.0
V = Ellipsoid(a,b,c)
Le = 1e-4
Lu = Inf #no magnetic diffusion
Ω = [0,0,1.0]
B₀ = [-y/b^2,x/a^2,0] #modified Malkus field.

#setup the problem to be solved
p = MHDProblem(N, V, Ω, Le, Lu, B₀, QGBasis, ConductingMFBasis)

#assemble the matrices using all threads available to julia (to start with 16 threads use "julia -t 16")
assemble!(p; threads=true)

#convert sparse to dense matrices and use LAPACK to solve generalized eigenproblem
esol = eigen(Matrix(p.RHS), Matrix(p.LHS))

#eigenvalues are complex with vanishingly small real part. Take only positive
#imaginary values:
evals = imag.(esol.values)[imag.(esol.values).>0]
evals = sort(evals,rev=true)

println("The ten largest eigen frequencies are:\n")
println(evals[1:10])

println("\n The ten smallest eigen frequencies are:\n")
println(reverse(evals[end-10:end]))
