var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Mire.jl Documentation",
    "title": "Mire.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#Mire.jl-Documentation-1",
    "page": "Mire.jl Documentation",
    "title": "Mire.jl Documentation",
    "category": "section",
    "text": "CurrentModule = Mire"
},

{
    "location": "#Mire.mat_force_galerkin!",
    "page": "Mire.jl Documentation",
    "title": "Mire.mat_force_galerkin!",
    "category": "function",
    "text": "mat_force_galerkin!(A,vs,N,forcefun,a,b,c, args...)\n\nFills Matrix A with Galerkin coefficients of force given by the function forcefun(u,a,b,c,args...).\n\n\n\n\n\n"
},

{
    "location": "#Functions-1",
    "page": "Mire.jl Documentation",
    "title": "Functions",
    "category": "section",
    "text": "mat_force_galerkin!"
},

{
    "location": "#Mire.assemblemhd",
    "page": "Mire.jl Documentation",
    "title": "Mire.assemblemhd",
    "category": "function",
    "text": "assemblemhd(N,a,b,c,Ω,b0)\n\nAssembles MHD eigen system, such that λAx=Bx\n\nThis is the dissipationless model, with\n\n∂ₜu = -2Ω×u + (∇×b0)×b + (∇×b)×b0 ∂ₜb = ∇×(u×b0)\n\nwith Ω = 1/Le * eΩ.\n\n\n\n\n\n"
},

{
    "location": "#Setting-up-the-eigen-problem-1",
    "page": "Mire.jl Documentation",
    "title": "Setting up the eigen problem",
    "category": "section",
    "text": "assemblemhd"
},

]}
