using Mire, Plots
revise()
a,b,c = 1.,1.,1.
Le = 1e-6
Ω = 1/Le*ez
n=4
B0 = [-y,x,0]
# vs = Mire.v(n,a,b,c);

LHS,RHS,vs = assemblemhd(n,a,b,c,Ω,B0)

esol = eigen(inv(Matrix(LHS))*RHS)

vtest = Mire.eigenvel(n,vs,esol.vectors,1,a,b,c)

plot(abs.(esol.values).+eps(),yscale=:log10,xscale=:log10)
