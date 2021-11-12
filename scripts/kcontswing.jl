using .NetworkResilience 
using PyPlot 
using Statistics 

# q=0, n=16, (1,15,0)
# q=0, n=32, (1,31,0)
# q=0, n=64, (1,63,0)

# q=0, n=16, (2,14,0)
# q=0, n=32, (2,30,0)
# q=0, n=64, (2,62,0)

# q=1, n=16, (1,15,0)
# q=1, n=32, (1,31,0)
# q=1, n=64, (1,63,0)

# q=0.1, n=16, (1,15,0)
# q=0.1, n=32, (1,31,0)
# q=0.1, n=64, (1,63,0)



n = 64
ns = 2
nd = 62
q = 1.0

κ = 5.0
I = 1.0
D = 1.0 

ensemblesize = 200
tol = 5e-4
acs = [] 
ks = []
fname = "data/kcontswing/q10_2_62_0.txt"
open(fname, "w") do f 
	for k=2:2:(n-3)
		println(k)
		ac, fmax, iter = swingcascadebisection(ensemblesize, tol, n, ns, nd, q, k, κ, I, D)
		push!(acs, ac) 
		push!(ks, k)
		x = ks[end]
		y = acs[end]
		write(f, "$x $y \n")  
	end
end










