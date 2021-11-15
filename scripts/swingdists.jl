using .NetworkResilience 

n = 100
ns = 50
nd = 50
q = 0.0
κ = 5.0
I = 1.0
D = 1.0 
tol = 5e-4
k = 4
ensemblesize = 1000
alpha = 2.0
acs = []

println("Running cascacde simulations for q=0")
for i=1:ensemblesize
	println(i) 
	ac, fmax, iter = swingcascadebisection(1, tol, n, ns, nd, q, k, κ, I, D)
	push!(acs, ac) 
end

fname = "data/hist_q0_50_50_0.txt"
open(fname, "w") do f 
	for i=1:length(acs)
		x = acs[i]
		write(f, "$x \n") 
	end
end

q = 0.1
acs = []

println("Running cascacde simulations for q=0.1")
for i=1:ensemblesize
	println(i) 
	ac, fmax, iter = swingcascadebisection(1, tol, n, ns, nd, q, k, κ, I, D)
	push!(acs, ac) 
end

fname = "data/hist_q1_50_50_0.txt"
open(fname, "w") do f 
	for i=1:length(acs)
		x = acs[i]
		write(f, "$x \n") 
	end
end