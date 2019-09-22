var documenterSearchIndex = {"docs":
[{"location":"man/examples/inertialmodes/#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"The examples are also available as notebooks in the source code folder example.","category":"page"},{"location":"man/examples/inertialmodes/#Inertial-modes-in-a-triaxial-ellipsoid-1","page":"Examples","title":"Inertial modes in a triaxial ellipsoid","text":"","category":"section"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"We want to solve the inertial mode equation","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"omega mathbfu = 2mathbfOmegatimesmathbfu","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"by expanding the velocity in a Cartesian polynomial basis and projecting onto these basis vectors following Lebovitz (1989).","category":"page"},{"location":"man/examples/inertialmodes/#Setting-up-the-problem-1","page":"Examples","title":"Setting up the problem","text":"","category":"section"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"The semi-axes of the triaxial ellipsoid","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"fracx^2a^2+fracy^2b^2+fracz^2c^2=1","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"are set to 1 for the sphere. Additionally we simplify the problem even further by taking the rotation axis along z, so that mathbfOmega=(001).","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"We truncate the problem at a maximum monomial x^iy^jz^k degree i+j+kleq N = 7.","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"# Include the packages in Julia\nusing Mire, LinearAlgebra, PyPlot #PyPlot uses matplotlib for plots","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"a,b,c = 1.1,0.9,0.7\nΩ = [0,0,1]\nN = 3","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"Assembling projects the left and right hand side of the inertial mode equation onto the basis vectors. For the integration a convenient formula is used (compare Lebovitz, 1989). Calling assemblehd outputs two sparse matrices A and B and the basis vectors uj. The Matrix B represents the left hand side and A the right hand side of","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"omega mathbfu = 2mathbfOmegatimesmathbfu","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"so that the eigen problem reads","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"omega Bmathbfu=Amathbfu","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"B,A, uj = assemblehd(N, a, b, c, Ω)","category":"page"},{"location":"man/examples/inertialmodes/#Solving-for-eigen-modes-1","page":"Examples","title":"Solving for eigen modes","text":"","category":"section"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"There are several ways to solve for eigen solutions of the generalized eigen problem. For small matrices we can simply invert Matrix B to reduce the problem to a standard eigen problem","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"B^-1Amathbfu=omegamathbfu","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"In Julia, the LAPACK routines for dense eigen problems are included in the standard library LinearAlgebra. Since A and B are sparse for now we have to convert B to a dense array by Matrix(B) before calling the inverse function inv. This is only feasible for small N, since we are now dealing with dense arrays. For larger N and thus larger matrices iterative sparse solvers should be applied.","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"esol = eigen(inv(Matrix(B))*A)","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"The eigen values and vectors are accessed by esol.values and esol.vectors respectively.","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"The eigen vectors contain the coefficients a_ji, so that the eigen velocity mathbfu_i is given by","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"mathbfu_i = sum_ja_jimathbfu_j","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"where mathbfu_j is the j-th basis vector in uj.","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"We can reconstruct the k-th eigenvelocity mathbfu_k by calling eigenvel:","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"k=length(esol.values)-3\nu_k = eigenvel(N,uj,esol.vectors,k,a,b,c)","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"u_k is now an array of cartesian polynomials with complex coefficients.","category":"page"},{"location":"man/examples/inertialmodes/#Plotting-the-mode-1","page":"Examples","title":"Plotting the mode","text":"","category":"section"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"include(joinpath(dirname(pathof(Mire)),\"../example/plotting.jl\"))\n\nfunction plotmode(a,b,c,u_k; kwargs...)\n    figure()\n    plot_velocity_equator(a,b,u_k; kwargs...)\n    title(\"x-y plane\")\n#     colorbar()\n    figure()\n    plot_velocity_meridional_x(b,c,u_k; kwargs...)\n    title(\"y-z plane\")\n#     colorbar()\n    figure()\n    plot_velocity_meridional_y(a,c,u_k; kwargs...)\n    title(\"x-z plane\")\n#     colorbar()\nend","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"We plot the m=2 quasi-geostrophic eigen mode with a frequency of","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"println(\"ω = \",imag.(esol.values[k]),\"𝕚\")","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"ω = 0.23780828249416838𝕚","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"plotmode(a,b,c,u_k, density=1.4, cmap=:plasma)","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"(Image: png)","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"(Image: png)","category":"page"},{"location":"man/examples/inertialmodes/#","page":"Examples","title":"Examples","text":"(Image: png)","category":"page"},{"location":"#Mire.jl-1","page":"Mire.jl","title":"Mire.jl","text":"","category":"section"},{"location":"#","page":"Mire.jl","title":"Mire.jl","text":"Modes In Rotating Ellipsoids","category":"page"},{"location":"#","page":"Mire.jl","title":"Mire.jl","text":"Toolbox written in Julia to solve eigen modes of rapidly rotating hydrodynamics or magnetohydrodynamics in a triaxial ellipsoid using Cartesian polynomials.","category":"page"},{"location":"#Development-1","page":"Mire.jl","title":"Development","text":"","category":"section"},{"location":"#","page":"Mire.jl","title":"Mire.jl","text":"Mire.jl is currently developed by Felix Gerick","category":"page"},{"location":"#Cite-1","page":"Mire.jl","title":"Cite","text":"","category":"section"},{"location":"#","page":"Mire.jl","title":"Mire.jl","text":"TODO","category":"page"},{"location":"man/functions/#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"CurrentModule = Mire","category":"page"},{"location":"man/functions/#Setting-up-the-eigen-problem-1","page":"Functions","title":"Setting up the eigen problem","text":"","category":"section"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"TODO: Description how to set up the problem here.","category":"page"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"assemblehd","category":"page"},{"location":"man/functions/#Mire.assemblehd","page":"Functions","title":"Mire.assemblehd","text":"assemblehd(N::Int, a::T, b::T, c::T, Ω ; dtype::DataType=BigFloat, kwargs...) where T\n\nAssemble the sparse matrices of the MHD mode problem. Returns right hand side `A`,\nleft hand side `B` and basis vectors `vs`.\n\n#Arguments:\n- `N`: maximum monomial degree\n- `a`: semi-axis x\n- `b`: semi-axis y\n- `c`: semi-axis z\n- `Ω`: rotation vector\n- `dtype`: datatype, default `BigFloat` for integration of monomials\n- `kwargs`: other keyword arguments passed to lower functions\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"assemblemhd\n","category":"page"},{"location":"man/functions/#Mire.assemblemhd","page":"Functions","title":"Mire.assemblemhd","text":"assemblemhd(N::Int, a::T, b::T, c::T, Ω, b0; dtype::DataType=BigFloat, kwargs...) where T\n\nAssemble the sparse matrices of the MHD mode problem. Returns right hand side A, left hand side B and basis vectors vs.\n\n#Arguments:\n\nN: maximum monomial degree\na: semi-axis x\nb: semi-axis y\nc: semi-axis z\nΩ: rotation vector\nb0: mean magnetic field vector\ndtype: datatype, default BigFloat for integration of monomials\nkwargs: other keyword arguments passed to lower functions\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#Low-level-functions-1","page":"Functions","title":"Low level functions","text":"","category":"section"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"Functions for low level control of the problem.","category":"page"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"projectforce","category":"page"},{"location":"man/functions/#Mire.projectforce","page":"Functions","title":"Mire.projectforce","text":"projectforce(N::Integer, vs, forcefun::Function, a::T, b::T, c::T, args...) where T\n\nAllocates new matrix A and fills elements by calling projectforce!(A,vs,N,forcefun,a,b,c, args...).\n\nCached version: projectforce(N,cmat,vs,forcefun,a,b,c, args...)\n\nwhere cmat[i,j,k] contains the integrals of monomials xⁱyʲzᵏ.\n\n#Arguments:\n\nN: maximum monomial degree\nvs: basis vectors\nforcefun: function of the force, e.g. coriolis\na: semi-axis x\nb: semi-axis y\nc: semi-axis z\nargs: other arguments needed for forcefun\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"projectforce!","category":"page"},{"location":"man/functions/#Mire.projectforce!","page":"Functions","title":"Mire.projectforce!","text":"projectforce!(A::AbstractArray{T, 2}, cmat::Array{T, 3}, vs::Array{Array{P, 1}, 1}, N::Integer, forcefun::Function, a::T, b::T, c::T, args...; kwargs...) where {T, P <: Polynomial{T}}\n\nDOCSTRING\n\n#Arguments:\n\nA: pre-allocated array\ncmat: pre-cached monomial integration values\nvs: basis vectors\nN: maximum monomial degree\nforcefun: function of the force, e.g. coriolis\na: semi-axis x\nb: semi-axis y\nc: semi-axis z\nargs: other arguments needed for forcefun\nkwargs: other keyword arguments\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"vel","category":"page"},{"location":"man/functions/#Mire.vel","page":"Functions","title":"Mire.vel","text":"vel(N,a,b,c)\n\nCompute all velocity basis vectors for a given maximal degree N (Lebovitz 1989, eq. 41-42)\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"eigenvel","category":"page"},{"location":"man/functions/#Mire.eigenvel","page":"Functions","title":"Mire.eigenvel","text":"eigenvel(N,vs,αs,a,b,c; norm=true)\n\nReconstructs velocity u, following Vidal & Cebron 2017 eq. (3.5).\n\nIf norm keyword is set true the velocity is normalised to satisfy int ucdot u dV=1.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"int_ellipsoid_surface","category":"page"},{"location":"man/functions/#Mire.int_ellipsoid_surface","page":"Functions","title":"Mire.int_ellipsoid_surface","text":"Integral over the surface of an ellipsoid.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"int_monomial_ellipsoid","category":"page"},{"location":"man/functions/#Mire.int_monomial_ellipsoid","page":"Functions","title":"Mire.int_monomial_ellipsoid","text":"int_monomial_ellipsoid(i::BigInt, j::BigInt, k::BigInt, a::T, b::T, c::T; dtype::DataType=BigFloat) where T\n\nDOCSTRING\n\n#Arguments:\n\ni: DESCRIPTION\nj: DESCRIPTION\nk: DESCRIPTION\na: DESCRIPTION\nb: DESCRIPTION\nc: DESCRIPTION\ndtype: DESCRIPTION\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"int_polynomial_ellipsoid","category":"page"},{"location":"man/functions/#Mire.int_polynomial_ellipsoid","page":"Functions","title":"Mire.int_polynomial_ellipsoid","text":"int_polynomial_ellipsoid(p, cmat)\n\nDOCSTRING\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#","page":"Functions","title":"Functions","text":"inner_product","category":"page"},{"location":"man/functions/#Mire.inner_product","page":"Functions","title":"Mire.inner_product","text":"inner_product(u,v,a,b,c)\n\nDefines inner product in an ellipsoidal volume intlangle uvrangle dV.\n\n\n\n\n\ninner_product(cmat, u, v; thresh=eps())\n\nDOCSTRING\n\n#Arguments:\n\ncmat: DESCRIPTION\nu: DESCRIPTION\nv: DESCRIPTION\nthresh: DESCRIPTION\n\n\n\n\n\n","category":"function"}]
}
