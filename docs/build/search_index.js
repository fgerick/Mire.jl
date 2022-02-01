var documenterSearchIndex = {"docs":
[{"location":"man/examples/modifiedmalkusmodes/#Modified-Malkus-modes-in-an-ellipsoid","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in an ellipsoid","text":"","category":"section"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"Some details are skipped here, which are already shown in the inertial mode example.","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"We want to solve the linearized momentum equation including the Lorentz force","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"partial_t mathbfu = -frac2mathrmLemathbfOmegatimesmathbfu-nabla p + left(nablatimesmathbfBright)timesmathbfB_0 + left(nablatimesmathbfB_0right)timesmathbfB","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"with mathbfB_0 = mathbfB_0(mathbfr), and leftmathbfumathbfBright = left mathbfumathbfBright(mathbfr)exp(iomega t). ","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"The evolution of the magnetic field is given by the diffusionless induction equation","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"partial_t mathbfB = nablatimesleft(mathbfutimesmathbfB_0right)","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"The magnetohydrodynamic problem is non-dimensionalized using a characteristic length L (e.g. core radius) and the Alfv\\'en wave period tau_A as a characteristic length scale. The arising non-dimensional number is the Lehnert number Le = B_0(LOmega sqrtrhomu_0), giving the ratio between the rotation period and the Alfv\\'en wave period.","category":"page"},{"location":"man/examples/modifiedmalkusmodes/#Setting-up-and-solving-the-problem","page":"Modified Malkus modes in Ellipsoid","title":"Setting up and solving the problem","text":"","category":"section"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"\nusing Mire, LinearAlgebra\n\na,b,c = 1.1,0.9,0.7\nV = Ellipsoid(a,b,c)\nΩ = [0.0,0.0,1.0]\nLe = 1e-4 #Lehnert number.\nN = 3\nB₀ = [-y/b^2,x/a^2,0] #Modified Malkus field adapted to the shape of the Ellipsoid.\n\n# create magnetohydrodynamic problem, using 3-D LebovitzBasis for the velocity \n# and ConductingMFBasis for the magnetic field:\np = MHDProblem(N, V, Ω, Le, B₀, LebovitzBasis, ConductingMFBasis) \nassemble!(p; threads=false)\n","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"Then solve the dense generalized eigen problem","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"A, B = Matrix(p.RHS), Matrix(p.LHS)\nevals, evecs = eigen(A, B)","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"Given the eigen values evals and eigen vectors evecs.","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"The eigen vectors mathbfx_i contain the n=n_u+n_b coefficients x_ji, so that the eigen velocity mathbfu_i and magnetic field mathbfB_i are given by","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"mathbfu_i = sum_j=1^n_ux_jitildemathbfu_j mathbfB_i = sum_j=n_u+1^nx_jitildemathbfB_j","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"where tildemathbfu_j, tildemathbfB_j are the j-th basis vectors in p.vbasis.el and p.bbasis.el.","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"We can reconstruct the velocities mathbfu_i and magnetic fields mathbfB_i for all i by calling velocities and magneticfields, respectively:","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"u = velocities(p.vbasis.el, evecs)\nB = magneticfields(p.bbasis.el, evecs)","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"u and B are now arrays of 3-D vectors where the components are Cartesian polynomials with complex coefficients.","category":"page"},{"location":"man/examples/modifiedmalkusmodes/#Plotting-the-kinetic-to-magnetic-energy-spectrum-of-modes","page":"Modified Malkus modes in Ellipsoid","title":"Plotting the kinetic to magnetic energy spectrum of modes","text":"","category":"section"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"We can calculate the kinetic energy ","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"E_mathrmkin = frac12intmathbfucdotmathbfumathrmdV","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"and magnetic energy","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"E_mathrmmag = frac12intmathbfBcdotmathbfBmathrmdV","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"for all modes by calling","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"ekin = [inner_product(u,u, p.cmat)/2 for u in u]\nemag = [inner_product(B,B, p.cmat)/2 for B in B]","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"Then, the frequency to Energy ratio spectrum can be plotted by","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"using PyPlot\n\nekm = abs.(ekin./emag)\n\nω = abs.(imag.(evals))\nloglog(ω, ekm, \".\")\nxlabel(L\"\\omega\\, [1/\\tau_A]\")\nylabel(L\"E_{kin}/E_{mag}\")","category":"page"},{"location":"man/examples/modifiedmalkusmodes/","page":"Modified Malkus modes in Ellipsoid","title":"Modified Malkus modes in Ellipsoid","text":"(Image: png)","category":"page"},{"location":"man/examples/inertialmodes/#Inertial-modes-in-an-ellipsoid","page":"Inertial modes in Ellipsoid","title":"Inertial modes in an ellipsoid","text":"","category":"section"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"We want to solve the inertial mode equation","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"partial_t mathbfu = -2mathbfOmegatimesmathbfu-nabla p","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"by expanding the velocity in a Cartesian polynomial basis and projecting onto these basis vectors following Lebovitz (1989).","category":"page"},{"location":"man/examples/inertialmodes/#Setting-up-the-problem","page":"Inertial modes in Ellipsoid","title":"Setting up the problem","text":"","category":"section"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"# Include the packages in Julia\nusing Mire, LinearAlgebra, PyPlot #PyPlot uses matplotlib for plots","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"The triaxial ellipsoid is defined by","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"fracx^2a^2+fracy^2b^2+fracz^2c^2=1","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"The rotation axis is taken to be along z, so that mathbfOmega=(001). We truncate at a maximum monomial polynomial degree N = 3, so that each monomial x^i y^j z^k has i+j+kleq N.","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"a,b,c = 1.1,0.9,0.7\nV = Ellipsoid(a,b,c) #volume\nΩ = [0.0,0.0,1.0] #rotation axis\nN = 3 #truncation degree\n\n# create hydrodynamic problem, using 3-D LebovitzBasis for the velocity:\np = HDProblem(N, V, Ω, LebovitzBasis) ","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"From here we can assemble the problem p. This means that we project the left and right hand side of the inertial mode equation onto the basis vectors mathbfu_j given by the LebovitzBasis. The pressure gradient force vanishes naturally in the projection, due to the incompressibility of the velocity. For the integration of the Cartesian polynomials (or rather the individual monomials) a convenient formula is used (compare Lebovitz, 1989). ","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"We can assemble the problem p by calling (threads=true enables multithreading to accelerate larger scale problems) ","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"assemble!(p; threads=false)","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"The left and right hand side matrices p.LHS and p.RHS, respectively, then represent the left and right hand side of","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"omega int mathbfu_i cdotmathbfu_j mathrmdV = -2int (mathbfOmegatimesmathbfu_i)cdotmathbfu_j mathrmdV","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"so that the eigen problem reads","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"omega Bmathbfx=Amathbfx","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"With B = p.LHS and A = p.RHS.","category":"page"},{"location":"man/examples/inertialmodes/#Solving-for-eigen-modes","page":"Inertial modes in Ellipsoid","title":"Solving for eigen modes","text":"","category":"section"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"There are several ways to solve for eigen solutions of the generalized eigen problem. For small matrices we can simply solve directly","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"omega Bmathbfx=Amathbfx","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"This can be done in Julia using the LAPACK routines for dense eigen problems are included in the standard library LinearAlgebra. Since A and B are sparse for now we have to convert B to a dense array by caling Matrix. This is only feasible for small N, since we are now dealing with dense arrays. For larger N and thus larger matrices iterative sparse solvers should be applied. An example for the sparse eigen solvers using Arpack.jl is given elsewhere.","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"A, B = Matrix(p.RHS), Matrix(p.LHS)\nevals, evecs = eigen(A, B)","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"Given the eigen values evals and eigen vectors evecs.","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"The eigen vectors mathbfx_i contain the coefficients x_ji, so that the eigen velocity mathbfu_i is given by","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"mathbfu_i = sum_jx_jitildemathbfu_j","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"where tildemathbfu_j is the j-th basis vector, given in p.vbasis.el[j].","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"We can reconstruct the eigenvelocities mathbfu_i for all i by calling velocities:","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"u = velocities(p.vbasis.el, evecs)","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"u is now an array of 3-D vectors where the components are Cartesian polynomials with complex coefficients.","category":"page"},{"location":"man/examples/inertialmodes/#Plotting-the-mode","page":"Inertial modes in Ellipsoid","title":"Plotting the mode","text":"","category":"section"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"An example of plotting streamlines at equatorial and meridional sections of one of the modes.","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"Include some PyPlot.jl plotting routines.","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"include(joinpath(dirname(pathof(Mire)),\"../example/plotting.jl\"))\n\nfunction plotmode(a,b,c,v_k; kwargs...)\n    figure()\n    plot_velocity_equator(a,b,v_k; kwargs...)\n    title(\"x-y plane\")\n#     colorbar()\n    figure()\n    plot_velocity_meridional_x(b,c,v_k; kwargs...)\n    title(\"y-z plane\")\n#     colorbar()\n    figure()\n    plot_velocity_meridional_y(a,c,v_k; kwargs...)\n    title(\"x-z plane\")\n#     colorbar()\nend","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"We plot the m=2 quasi-geostrophic eigen mode with a frequency of","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"k = findfirst(0.23 .< abs.(imag.(evals)) .< 0.24)\nprintln(\"ω = \",imag.(evals[k]),\"𝕚\")\n#ω = -0.23780828249417196𝕚","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"And to plot ","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"plotmode(a,b,c,u[k], density=1.4, cmap=:plasma)","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"(Image: png)","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"(Image: png)","category":"page"},{"location":"man/examples/inertialmodes/","page":"Inertial modes in Ellipsoid","title":"Inertial modes in Ellipsoid","text":"(Image: png)","category":"page"},{"location":"#Mire.jl","page":"Home","title":"Mire.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modes In Rotating Ellipsoids","category":"page"},{"location":"","page":"Home","title":"Home","text":"Toolbox written in Julia to solve eigen modes of rapidly rotating (magneto)hydrodynamics in an ellipsoid using Cartesian polynomials.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Some example usage is given in the section Documentation and in the Continuous Integration file.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Install Mire.jl from the julia REPL prompt with","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add https://github.com/fgerick/Mire.jl.git","category":"page"},{"location":"#Development","page":"Home","title":"Development","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Mire.jl is currently developed by Felix Gerick","category":"page"},{"location":"man/functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"CurrentModule = Mire","category":"page"},{"location":"man/functions/#Defining-the-Volume","page":"Functions","title":"Defining the Volume","text":"","category":"section"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"Ellipsoid","category":"page"},{"location":"man/functions/#Mire.Ellipsoid","page":"Functions","title":"Mire.Ellipsoid","text":"Ellipsoid{T<:Number} <: Volume{T}\n\nVolume type Ellipsoid. Create by calling Ellipsoid(a,b,c), with a,b,c the semi-axes. a,b,c can be of any Number type.\n\nExamples: Ellipsoid(1.1,1.0,0.9) is a Ellipsoid{Float64}, Ellipsoid(1//1,1//2,1//5) is a Ellispoid{Rational{Int64}}.\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"Sphere","category":"page"},{"location":"man/functions/#Mire.Sphere","page":"Functions","title":"Mire.Sphere","text":"Sphere{T<:Number} <: Volume{T}\n\nVolume type Sphere. Create by calling Sphere{T}(), with T any Number type. Default: Sphere(T) gives a Sphere{Float64}(). For other types use, e.g. Sphere{Rational{BigInt}}().\n\n\n\n\n\n","category":"type"},{"location":"man/functions/#Vector-bases","page":"Functions","title":"Vector bases","text":"","category":"section"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"LebovitzBasis","category":"page"},{"location":"man/functions/#Mire.LebovitzBasis","page":"Functions","title":"Mire.LebovitzBasis","text":"LebovitzBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}\n\nBasis of 3-D vector field, so that mathbfucdotmathbfn = 0 at partialmathcalV and nablacdotmathbfu = 0 after Lebovitz (1989).\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"QGBasis","category":"page"},{"location":"man/functions/#Mire.QGBasis","page":"Functions","title":"Mire.QGBasis","text":"QGBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}\n\nBasis of QG vector field (Gerick et al., 2020), so that mathbfu=nabla(h^3x^ny^m)timesnabla(zh).\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"QGIMBasis","category":"page"},{"location":"man/functions/#Mire.QGIMBasis","page":"Functions","title":"Mire.QGIMBasis","text":"QGIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}\n\nBasis of complex QG inertial modes (Maffei et al., 2017). This basis is orthonormal.\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"QGRIMBasis","category":"page"},{"location":"man/functions/#Mire.QGRIMBasis","page":"Functions","title":"Mire.QGRIMBasis","text":"QGRIMBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}\n\nReal basis similar to the QG inertial modes (Maffei et al., 2017). They are not the solutions to the QG inertial mode equation and they are not orthogonal.\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"ConductingMFBasis","category":"page"},{"location":"man/functions/#Mire.ConductingMFBasis","page":"Functions","title":"Mire.ConductingMFBasis","text":"ConductingMFBasis{T<:Number,Vol<:Volume{T}} <: VectorBasis{T,Vol}\n\nBasis of 3-D magnetic field, so that mathbfBcdotmathbfn = 0 at partialmathcalV and nablacdotmathbfB = 0 after Lebovitz (1989). It is exactly the same as LebovitzBasis.\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"InsulatingMFBasis","category":"page"},{"location":"man/functions/#Mire.InsulatingMFBasis","page":"Functions","title":"Mire.InsulatingMFBasis","text":"InsulatingMFBasis{T<:Number,Vol<:Sphere{T}} <: VectorBasis{T,Vol}\n\nBasis of insulating magnetic fields following Gerick et al. (2021). For now, only for Vol<:Sphere{T}, i.e. in a spherical domain. \n\n\n\n\n\n","category":"type"},{"location":"man/functions/#Setting-up-the-problem","page":"Functions","title":"Setting up the problem","text":"","category":"section"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"HDProblem","category":"page"},{"location":"man/functions/#Mire.HDProblem","page":"Functions","title":"Mire.HDProblem","text":"HDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T, Vol}\n\nDefines hydrodynamic problem.\n\nExample:\n\nN = 5\nΩ = [0,0,1.0]\nV = Ellipsoid(1.1,1.0,0.9)\nproblem = HDProblem(N,V,Ω,LebovitzBasis)\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"MHDProblem","category":"page"},{"location":"man/functions/#Mire.MHDProblem","page":"Functions","title":"Mire.MHDProblem","text":"MHDProblem{T<:Number,Vol<:Volume{T}} <: MireProblem{T,Vol}\n\nDefines magnetohydrodynamic problem.\n\nExample for a hybrid QG model with 3-D magnetic field with  conducting boundary condition and QG velocities:\n\nN = 5\nΩ = [0,0,1.0]\na,b,c = 1.1,1.0,0.9\nV = Ellipsoid(a,b,c)\nB0 = [-y/b^2,x/a^2,0] #Malkus field\nproblem = MHDProblem(N,V,Ω,B0,QGBasis,ConductingMFBasis)\n\n\n\n\n\n","category":"type"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"assemble!","category":"page"},{"location":"man/functions/#Mire.assemble!","page":"Functions","title":"Mire.assemble!","text":"assemble!(P::HDProblem{T,V}; threads=false, verbose=false, kwargs...) where {T,V}\n\nAssembles the matrices P.LHS and P.RHS, i.e. projecting the velocity basis P.vbasis on the inertial acceleration and Coriolis force.\n\n\n\n\n\nassemble!(P::MHDProblem{T,V}) where {T,V}\n\nAssembles the matrices P.LHS and P.RHS, i.e. projecting the velocity and magnetic field bases on the inertial acceleration, Coriolis force, Lorentz force and mgnetic advection.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/#Low-level-functions","page":"Functions","title":"Low level functions","text":"","category":"section"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"Functions for low level control of the problem.","category":"page"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"projectforce","category":"page"},{"location":"man/functions/#Mire.projectforce","page":"Functions","title":"Mire.projectforce","text":"projectforce(vs_i, vs_j, cmat, forcefun, args...; kwargs...)\n\nProject basis vs_i onto the forcing forcefun(vs_j, args...) using precached monomials in cmat.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"projectforce!","category":"page"},{"location":"man/functions/#Mire.projectforce!","page":"Functions","title":"Mire.projectforce!","text":"projectforce!(i0, j0, itemps, jtemps, valtemps, cmat, vs_i, vs_j, forcefun, args...; verbose=false, thresh=10eps())\n\nProject basis vs_i onto the forcing forcefun(vs_j, args...). Pushes indices and values of non-zero entries into the itemps,jtemps and valtems arrays. i0 and j0 are the origin of indices (useful for combining multiple equations, should be 0 to have no shift in indices). See projectforce to construct a sparse array from the indices and values.\n\n\n\n\n\n","category":"function"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"projectforcet!","category":"page"},{"location":"man/functions/#Mire.projectforcet!","page":"Functions","title":"Mire.projectforcet!","text":"projectforcet!(args...)\n\nMultithreaded version of projectforce!. Input has to be adapted to number of threads!\n\n\n\n\n\n","category":"function"},{"location":"man/functions/","page":"Functions","title":"Functions","text":"inner_product","category":"page"},{"location":"man/functions/#Mire.inner_product","page":"Functions","title":"Mire.inner_product","text":"inner_product(u, v, cmat)\n\nInner product in an ellipsoidal volume int  mathbfucdot mathbfv dV,  using precached monomial integrations in cmat.\n\n\n\n\n\ninner_product(u, v, a, b, c)\n\nInner product in an ellipsoidal volume int  mathbfucdot mathbfv dV.\n\nwarning: Missing factor\nThe factor pi is removed in the integration to be able to use exact integration using Rational.  Remember to reintroduce it when the actual value of the integration is needed!\n\n\n\n\n\n","category":"function"}]
}
